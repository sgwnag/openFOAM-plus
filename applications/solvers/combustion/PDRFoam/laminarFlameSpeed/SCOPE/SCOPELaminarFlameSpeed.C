/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "IFstream.H"
#include "SCOPELaminarFlameSpeed.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarFlameSpeedModels
{
    defineTypeNameAndDebug(SCOPE, 0);

    addToRunTimeSelectionTable
    (
        laminarFlameSpeed,
        SCOPE,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

static void readPolynomialCoeffs
(
    Polynomial<7>& coeffs,
    const word& keyword,
    const dictionary& dict
)
{
    {
        coeffs[0] = dict.get<scalar>(keyword + "0");
        coeffs[1] = dict.get<scalar>(keyword + "1");
        coeffs[2] = dict.get<scalar>(keyword + "2");
        coeffs[3] = dict.get<scalar>(keyword + "3");
        coeffs[4] = dict.get<scalar>(keyword + "4");
        coeffs[5] = dict.get<scalar>(keyword + "5");
        coeffs[6] = dict.get<scalar>(keyword + "6");
    }

    // TBD: support direct reading of all coeffs?
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::SCOPE::polynomial::polynomial
(
    const dictionary& polyDict
)
:
    FixedList<scalar, 7>(polyDict.lookup("coefficients")),
    ll(polyDict.get<scalar>("lowerLimit")),
    ul(polyDict.get<scalar>("upperLimit")),
    llv(polyPhi(ll, *this)),
    ulv(polyPhi(ul, *this)),
    lu(0)
{}


Foam::laminarFlameSpeedModels::SCOPE::SCOPE
(
    const dictionary& dict,
    const psiuReactionThermo& ct
)
:
    laminarFlameSpeed(dict, ct),

    coeffsDict_
    (
        dictionary
        (
            IFstream(dict.get<fileName>("fuelFile"))()
        ).subDict(typeName + "Coeffs")
    ),
    LFL_(coeffsDict_.get<scalar>("lowerFlamabilityLimit")),
    UFL_(coeffsDict_.get<scalar>("upperFlamabilityLimit")),
    SuPolyL_(coeffsDict_.subDict("lowerSuPolynomial")),
    SuPolyU_(coeffsDict_.subDict("upperSuPolynomial")),
    Texp_(),
    pexp_(),
    CIn_(coeffsDict_.getOrDefault<scalar>("CIn", 0)),
    MaPolyL_(coeffsDict_.subDict("lowerMaPolynomial")),
    MaPolyU_(coeffsDict_.subDict("upperMaPolynomial"))
{
    readPolynomialCoeffs(Texp_, "Texp", coeffsDict_);
    readPolynomialCoeffs(pexp_, "pexp", coeffsDict_);

    SuPolyL_.ll = max(SuPolyL_.ll, LFL_) + SMALL;
    SuPolyU_.ul = min(SuPolyU_.ul, UFL_) - SMALL;

    SuPolyL_.lu = 0.5*(SuPolyL_.ul + SuPolyU_.ll);
    SuPolyU_.lu = SuPolyL_.lu - SMALL;

    MaPolyL_.lu = 0.5*(MaPolyL_.ul + MaPolyU_.ll);
    MaPolyU_.lu = MaPolyL_.lu - SMALL;

    if (debug)
    {
        Info<< "phi     Su  (T = Tref, p = pref)" << endl;
        const label n = 200;
        for (int i=0; i<n; i++)
        {
            scalar phi = (2.0*i)/n;
            Info<< phi << token::TAB << SuRef(phi) << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laminarFlameSpeedModels::SCOPE::~SCOPE()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::polyPhi
(
    scalar phi,
    const polynomial& a
)
{
    const scalar x = phi - 1.0;

    return
        a[0]
       *(
           scalar(1)
         + x*(a[1] + x*(a[2] + x*(a[3] + x*(a[4] + x*(a[5] + x*a[6])))))
        );
}


inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::SuRef
(
    scalar phi
) const
{
    if (phi < LFL_ || phi > UFL_)
    {
        // Return 0 beyond the flammability limits
        return scalar(0);
    }
    else if (phi < SuPolyL_.ll)
    {
        // Use linear interpolation between the low end of the
        // lower polynomial and the lower flammability limit
        return SuPolyL_.llv*(phi - LFL_)/(SuPolyL_.ll - LFL_);
    }
    else if (phi > SuPolyU_.ul)
    {
        // Use linear interpolation between the upper end of the
        // upper polynomial and the upper flammability limit
        return SuPolyU_.ulv*(UFL_ - phi)/(UFL_ - SuPolyU_.ul);
    }
    else if (phi < SuPolyL_.lu)
    {
        // Evaluate the lower polynomial
        return polyPhi(phi, SuPolyL_);
    }
    else if (phi > SuPolyU_.lu)
    {
        // Evaluate the upper polynomial
        return polyPhi(phi, SuPolyU_);
    }
    else
    {
        FatalErrorInFunction
            << "phi = " << phi
            << " cannot be handled by SCOPE function with the "
               "given coefficients"
            << exit(FatalError);

        return scalar(0);
    }
}

inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::Ma
(
    scalar phi
) const
{
    if (phi < LFL_ || phi > UFL_)
    {
        // Return 0 beyond the flamibility limits
        return scalar(0);
    }
    else if (phi < MaPolyL_.ll)
    {
        // Use linear interpolation between the low end of the
        // lower polynomial and the lower flammability limit
        return MaPolyL_.llv*(phi - LFL_)/(MaPolyL_.ll - LFL_);
    }
    else if (phi > MaPolyU_.ul)
    {
        // Use linear interpolation between the upper end of the
        // upper polynomial and the upper flammability limit
        return MaPolyU_.ulv*(UFL_ - phi)/(UFL_ - MaPolyU_.ul);
    }
    else if (phi < MaPolyL_.lu)
    {
        // Evaluate the lower polynomial
        return polyPhi(phi, MaPolyL_);
    }
    else if (phi > MaPolyU_.lu)
    {
        // Evaluate the upper polynomial
        return polyPhi(phi, MaPolyU_);
    }
    else
    {
        FatalErrorInFunction
            << "phi = " << phi
            << " cannot be handled by SCOPE function with the "
               "given coefficients"
            << exit(FatalError);

        return scalar(0);
    }
}


inline Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    scalar p,
    scalar Tu,
    scalar phi
) const
{
    constexpr scalar Tref = 300.0;
    constexpr scalar pRef = 1.013e5;

    const scalar Texp = Texp_.value(phi-1.0);
    const scalar pexp = pexp_.value(phi-1.0);

    return SuRef(phi)*pow((Tu/Tref), Texp)*pow((p/pRef), pexp);
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    scalar phi
) const
{
    auto tSu0 = tmp<volScalarField>::New
    (
        IOobject
        (
            "Su0",
            p.time().timeName(),
            p.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p.mesh(),
        dimensionedScalar(dimVelocity, Zero)
    );
    auto& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi);
    }

    forAll(Su0.boundaryField(), patchi)
    {
        scalarField& Su0p = Su0.boundaryFieldRef()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tup = Tu.boundaryField()[patchi];

        forAll(Su0p, facei)
        {
            Su0p[facei] = Su0pTphi(pp[facei], Tup[facei], phi);
        }
    }

    return tSu0;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Su0pTphi
(
    const volScalarField& p,
    const volScalarField& Tu,
    const volScalarField& phi
) const
{
    auto tSu0 = tmp<volScalarField>::New
    (
        IOobject
        (
            "Su0",
            p.time().timeName(),
            p.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        p.mesh(),
        dimensionedScalar(dimVelocity, Zero)
    );
    auto& Su0 = tSu0.ref();

    forAll(Su0, celli)
    {
        Su0[celli] = Su0pTphi(p[celli], Tu[celli], phi[celli]);
    }

    forAll(Su0.boundaryField(), patchi)
    {
        scalarField& Su0p = Su0.boundaryFieldRef()[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& Tup = Tu.boundaryField()[patchi];
        const scalarField& phip = phi.boundaryField()[patchi];

        forAll(Su0p, facei)
        {
            Su0p[facei] =
                Su0pTphi
                (
                    pp[facei],
                    Tup[facei],
                    phip[facei]
                );
        }
    }

    return tSu0;
}


Foam::scalar Foam::laminarFlameSpeedModels::SCOPE::CIn() const noexcept
{
    return CIn_ ;
}


Foam::tmp<Foam::volScalarField> Foam::laminarFlameSpeedModels::SCOPE::Ma
(
    const volScalarField& phi
) const
{
    auto tMa = tmp<volScalarField>::New
    (
        IOobject
        (
            "Ma",
            phi.time().timeName(),
            phi.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi.mesh(),
        dimensionedScalar(dimless, Zero)
    );
    auto& ma = tMa.ref();

    forAll(ma, celli)
    {
        ma[celli] = Ma(phi[celli]);
    }

    forAll(ma.boundaryField(), patchi)
    {
        scalarField& map = ma.boundaryFieldRef()[patchi];
        const scalarField& phip = phi.boundaryField()[patchi];

        forAll(map, facei)
        {
            map[facei] = Ma(phip[facei]);
        }
    }

    return tMa;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::SCOPE::Ma() const
{
    if (psiuReactionThermo_.composition().contains("ft"))
    {
        const volScalarField& ft = psiuReactionThermo_.composition().Y("ft");

        return Ma
        (
            dimensionedScalar
            (
                "stoichiometricAirFuelMassRatio",
                psiuReactionThermo_
            )*ft/(scalar(1) - ft)
        );
    }
    else
    {
        const fvMesh& mesh = psiuReactionThermo_.p().mesh();

        return tmp<volScalarField>::New
        (
            IOobject
            (
                "Ma",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("Ma", dimless, Ma(equivalenceRatio_))
        );
    }
}


Foam::tmp<Foam::volScalarField>
Foam::laminarFlameSpeedModels::SCOPE::operator()() const
{
    if (psiuReactionThermo_.composition().contains("ft"))
    {
        const volScalarField& ft = psiuReactionThermo_.composition().Y("ft");

        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            dimensionedScalar
            (
                "stoichiometricAirFuelMassRatio",
                psiuReactionThermo_
            )*ft/(scalar(1) - ft)
        );
    }
    else
    {
        return Su0pTphi
        (
            psiuReactionThermo_.p(),
            psiuReactionThermo_.Tu(),
            equivalenceRatio_
        );
    }
}


// ************************************************************************* //

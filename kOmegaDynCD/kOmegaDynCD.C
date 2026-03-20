/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "kOmegaDynCD.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void kOmegaDynCD<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = (k_/(omega_+low_omega_))*(cBeta_*Cmu_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDynCD<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime
        )
    );
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDynCD<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kOmegaDynCD<BasicMomentumTransportModel>::kOmegaDynCD
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.072
        )
    ),

    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),

    gamma_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.52
        )
    ),

    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            0.5
        )
    ),

    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),

    cBetaMin_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cBetaMin",
            this->coeffDict_,
            0.0
        )
    ),

    cBetaMax_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cBetaMax",
            this->coeffDict_,
            22.22
        )
    ),

    cBetaLim_
    (
        "cBetaLim",
        dimless,
        1.0/betaStar_.value()
    ),

    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",         
            this->runTime_.constant(),    
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    finestra_(
        this->coeffDict().lookupOrDefault("nPrecedents", 2)
    ),

    low_omega_
    (
        "low_omega",
        dimensionSet(0, 0, -1, 0, 0),
        SMALL
    ),

    low_k_
    (
        "low_k",
        dimensionSet(0, 2, -2, 0, 0),
        SMALL
    ),

    low_Z_
    (
        "low_Z",
        dimensionSet(0, 4, -4, 0, 0),
        SMALL
    ),

    low_cBeta_
    (
        "low_cBeta",
        dimless,
        1e-16
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    cBeta_
    (
        IOobject
        (
            "cBeta",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    z_
    (
        IOobject
        (
            "z",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("z", dimless, 1.0)
    ),

    S 
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedSymmTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    ),

    gU 
    (
        IOobject
        (
            "gU",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero)
    )

{
    Info << "kOmegaDynCD constructor called!" << nl;

////////////////////////////////    PtrLists    ////////////////////////////////////
    
    kPrev_.setSize(finestra_ - 1);
    omegaPrev_.setSize(finestra_ - 1);
    UPrev_.setSize(finestra_ - 1);
    kprint_.setSize(finestra_ - 1);
    omegaprint_.setSize(finestra_ - 1);
    Uprint_.setSize(finestra_ - 1);
    gUprint_.setSize(finestra_ - 1);

    for (int i = 0; i < finestra_ - 1; ++i)
    {
        kPrev_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
            )
        );

        omegaPrev_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "omegaPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,0,-1,0,0), SMALL)
            )
        );

        UPrev_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "UPrev_" + name(i),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
            )
        );
        kprint_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "kprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,2,-2,0,0), 0)
            )
        );

        omegaprint_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "omegaprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedScalar(dimensionSet(0,0,-1,0,0), 1e-12)
            )
        );

        Uprint_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "Uprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedVector(dimensionSet(0,1,-1,0,0), Zero)
            )
        );

        gUprint_.set
        (
            i,
            new volTensorField
            (
                IOobject
                (
                    "gUprint_" + name(i+1),
                    this->runTime_.timeName(),
                    this->mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                this->mesh_,
                dimensionedTensor(dimensionSet(0,0,-1,0,0), Zero)
            )
        );
    }

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kOmegaDynCD<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        betaStar_.readIfPresent(this->coeffDict());

        cBetaMin_.readIfPresent(this->coeffDict());
        cBetaMax_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicMomentumTransportModel>
void kOmegaDynCD<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local Ref. 
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;

    volScalarField& nut = this->nut_;

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));

    // Update omega and G at walls
    omega_.boundaryFieldRef().updateCoeffs();   

    // Specific dissipation equation
    tmp<fvScalarMatrix> omegaEqn(
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*G*omega_/(k_ + low_k_)
      - fvm::SuSp((2.0/3.0)*gamma_*alpha*rho*divU, omega_)
      - fvm::Sp(z_*(5.0/6.0)*(1/(cBeta_+low_cBeta_))*alpha*rho*omega_, omega_)
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(omega_);
    bound(omega_, this->omegaMin_);

    Info << "min(omega_): " << gMin(omega_) 
     << "  max(omega_): " << gMax(omega_) << nl;
    Info << "min(k_): " << gMin(k_)
     << "  max(k_): " << gMax(k_) << nl;

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn(
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(z_*(1/(cBeta_+low_cBeta_))*alpha*rho*omega_, k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

//////////////////////////////////  Dynamic procedure  ////////////////////////////////////

    const volTensorField& gradU = tgradU();

    S = symm(gradU);
    
    gU = gradU;

    volSymmTensorField sum_uu = symm(U * U);
    volVectorField sum_u = U;

    volSymmTensorField kES_mean = (2*Cmu_) * (k_/(omega_ + low_omega_)) * S;
    volSymmTensorField sum_S = S;

    volScalarField k_mean = k_;
    volScalarField omega_mean = omega_;

    volScalarField squaredS_mean = sqr(tr(S)) - sqr(tr(S)); // init placeholder (si potrebbe scrivere meglio)
    
    forAll(S, cellI)
    {
          squaredS_mean[cellI] =
          sqr(S[cellI].xx())
        + sqr(S[cellI].yy())
        + sqr(S[cellI].zz());
    }

    tgradU.clear();

    for (int i = 0; i < finestra_ - 1; ++i)
    {
        sum_uu += symm(UPrev_[i] * UPrev_[i]);
        sum_u += UPrev_[i];

        volTensorField gradUPrev_ = fvc::grad(UPrev_[i]);
        volSymmTensorField SPrev_ = symm(gradUPrev_);

        kES_mean += (2*Cmu_)*(kPrev_[i]/(omegaPrev_[i] + low_omega_)) * SPrev_;
        sum_S += SPrev_;

        k_mean += kPrev_[i];
        omega_mean += omegaPrev_[i];

        forAll(SPrev_, cellI)
        {
          squaredS_mean[cellI] +=
          sqr(SPrev_[cellI].xx())
		  + sqr(SPrev_[cellI].yy())
		  + sqr(SPrev_[cellI].zz());
        }
    }

    sum_u = sum_u / finestra_;
    volSymmTensorField uu_average = sum_uu / finestra_;
    volSymmTensorField u_average_u_average = symm((sum_u) * (sum_u));

    volSymmTensorField L = uu_average - u_average_u_average;

    kES_mean = kES_mean / finestra_;
    volSymmTensorField S_average = sum_S / finestra_;

    k_mean = k_mean / finestra_;
    volScalarField k_mean_T = k_mean + 0.5*(tr(uu_average) - (sum_u & sum_u)); 

    omega_mean = omega_mean / finestra_;

    squaredS_mean = squaredS_mean / finestra_;

    volScalarField SSquared_mean = sqr(tr(S_average)) - sqr(tr(S_average));
    
    forAll(S_average, cellI)
    {
        SSquared_mean[cellI] +=
            sqr(S_average[cellI].xx())
          + sqr(S_average[cellI].yy())
          + sqr(S_average[cellI].zz());
    }

    omega_mean = (1 / (k_mean_T + low_k_)) * ((k_mean * omega_mean) + ((this->nu()) * cBeta_) * (squaredS_mean - SSquared_mean));

    volSymmTensorField Z = kES_mean - (2*Cmu_) * (k_mean_T / (omega_mean + low_omega_)) * S_average;

    cBeta_ = (Z && L) / ((Z && Z) + low_Z_);

    //Debug 
    Info << "max value of cBeta: " << max(cBeta_) << nl;
	Info << "min value of cBeta: " << min(cBeta_) << nl;

    if (this->runTime_.timeIndex() < 10)
    {
        cBeta_ = 11.11;
    }

    // Clipping cBeta_ using user-settable limits
    const scalar CmuMin = cBetaMin_.value();
    const scalar CmuMax = cBetaMax_.value();

    forAll(cBeta_, cellI)
    {
        if (cBeta_[cellI] > CmuMax)
        {
            cBeta_[cellI] = CmuMax;
        }
        else if (cBeta_[cellI] < CmuMin)
        {
            cBeta_[cellI] = CmuMin;
        }
    }

    z_ = cBeta_ / cBetaLim_;

    // Update stored fields (previous-step)

    for (int i = 0; i <= finestra_ - 2; ++i)
    {
        Uprint_[i] = UPrev_[i];
        kprint_[i] = kPrev_[i];
        omegaprint_[i] = omegaPrev_[i];
        gUprint_[i] = fvc::grad(UPrev_[i]);
    }

    // shift history and insert current values at the end

    for (int i = 0; i < finestra_ - 2; ++i)
    {
        bound(kPrev_[i], this->kMin_);
        bound(omegaPrev_[i], this->omegaMin_);

        UPrev_[i] = UPrev_[i+1];
        kPrev_[i] = kPrev_[i+1];
        omegaPrev_[i] = omegaPrev_[i+1];
    }

    UPrev_[finestra_ - 2] = this->U_;
    kPrev_[finestra_ - 2] = this->k_;
    omegaPrev_[finestra_ - 2] = this->omega_;

    correctNut();

    Info << "Time = " << this->runTime_.timeName() << nl << endl;
}

////////////////////////////////////////////////////////////////////////////////

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
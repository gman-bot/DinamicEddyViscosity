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

#include "kOmegaDyn.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
void kOmegaDyn<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = (k_/(omega_+low_omega_))*(dynamicCmu_/betaStar_);
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDyn<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix(k_, dimVolume*this->rho_.dimensions()*k_.dimensions()/dimTime)
    );
}

template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kOmegaDyn<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>(
        new fvScalarMatrix(omega_, dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime)
    );
}

template<class BasicMomentumTransportModel>
kOmegaDyn<BasicMomentumTransportModel>::kOmegaDyn
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
    eddyViscosity<BasicMomentumTransportModel>(type, alpha, rho, U, alphaRhoPhi, phi, viscosity),
    beta_("beta", this->coeffDict_, 0.072),
    betaStar_("betaStar", this->coeffDict_, 0.09),
    gamma_("gamma", this->coeffDict_, 0.52),
    alphaK_("alphaK", this->coeffDict_, 0.5),
    alphaOmega_("alphaOmega", this->coeffDict_, 0.5),
    dynamicCmuMin_("dynamicCmuMin", this->coeffDict_, 0.0),
    dynamicCmuMax_("dynamicCmuMax", this->coeffDict_, 0.2),
    physicalProperties_(IOobject("physicalProperties", this->runTime_.constant(), this->mesh_, IOobject::MUST_READ, IOobject::NO_WRITE)),
    finestra_(this->coeffDict().lookupOrDefault("nPrecedents", 2)),
    low_omega_("low_omega", dimensionSet(0, 0, -1, 0, 0), SMALL),
    low_k_("low_k", dimensionSet(0, 2, -2, 0, 0), SMALL),
    low_Z_("low_Z", dimensionSet(0, 4, -4, 0, 0), SMALL),
    k_(IOobject(IOobject::groupName("k", alphaRhoPhi.group()), this->runTime_.timeName(), this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE), this->mesh_),
    omega_(IOobject(IOobject::groupName("omega", alphaRhoPhi.group()), this->runTime_.timeName(), this->mesh_, IOobject::MUST_READ, IOobject::AUTO_WRITE), this->mesh_),
    dynamicCmu_(IOobject("dynamicCmu", this->runTime_.timeName(), this->mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE), this->mesh_),
    S(IOobject("S", this->runTime_.timeName(), this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE), this->mesh_, dimensionedSymmTensor(dimensionSet(0, 0, -1, 0, 0), Zero)),
    gU(IOobject("gU", this->runTime_.timeName(), this->mesh_, IOobject::NO_READ, IOobject::NO_WRITE), this->mesh_, dimensionedTensor(dimensionSet(0, 0, -1, 0, 0), Zero))
{
    const int nPrev = finestra_ - 1;
    kPrev_.setSize(nPrev);
    omegaPrev_.setSize(nPrev);
    UPrev_.setSize(nPrev);
    SPrev_.setSize(nPrev);
    
    kprint_.setSize(nPrev);
    omegaprint_.setSize(nPrev);
    Uprint_.setSize(nPrev);
    gUprint_.setSize(nPrev);

    for (int i = 0; i < nPrev; ++i)
    {
        kPrev_.set(i, new volScalarField(IOobject("kPrev_"+name(i), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedScalar(dimViscosity, 0)));
        omegaPrev_.set(i, new volScalarField(IOobject("omegaPrev_"+name(i), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedScalar(low_omega_.dimensions(), SMALL)));
        UPrev_.set(i, new volVectorField(IOobject("UPrev_"+name(i), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedVector(dimVelocity, Zero)));
        SPrev_.set(i, new volSymmTensorField(IOobject("SPrev_"+name(i), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedSymmTensor(low_omega_.dimensions(), Zero)));
        
        kprint_.set(i, new volScalarField(IOobject("kprint_"+name(i+1), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedScalar(dimViscosity, 0)));
        omegaprint_.set(i, new volScalarField(IOobject("omegaprint_"+name(i+1), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedScalar(low_omega_.dimensions(), 1e-12)));
        Uprint_.set(i, new volVectorField(IOobject("Uprint_"+name(i+1), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedVector(dimVelocity, Zero)));
        gUprint_.set(i, new volTensorField(IOobject("gUprint_"+name(i+1), this->runTime_.timeName(), this->mesh_), this->mesh_, dimensionedTensor(low_omega_.dimensions(), Zero)));
    }

    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);
}

template<class BasicMomentumTransportModel>
bool kOmegaDyn<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<BasicMomentumTransportModel>::read())
    {
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
        dynamicCmuMin_.readIfPresent(this->coeffDict());
        dynamicCmuMax_.readIfPresent(this->coeffDict());
        return true;
    }
    return false;
}

template<class BasicMomentumTransportModel>
void kOmegaDyn<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_) return;

    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;

    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(this->mesh_));

    eddyViscosity<BasicMomentumTransportModel>::correct();

    volScalarField divU(fvc::div(fvc::absolute(this->phi(), U)));
    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField G(this->GName(), nut*(tgradU() && dev(twoSymm(tgradU()))));

    omega_.boundaryFieldRef().updateCoeffs();   

    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha*rho*G*omega_/(k_ + low_k_)
      - fvm::SuSp((2.0/3.0)*gamma_*alpha*rho*divU, omega_)
      - fvm::Sp(beta_*alpha*rho*omega_, omega_)
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    solve(omegaEqn);
    bound(omega_, this->omegaMin_);

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha*rho*divU, k_)
      - fvm::Sp(betaStar_*alpha*rho*omega_, k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    bound(k_, this->kMin_);

    // Procedura Dinamica Ottimizzata
    S = symm(tgradU());
    gU = tgradU();
    tgradU.clear();

    volSymmTensorField sum_uu = symm(U * U);
    volVectorField sum_u = U;
    volSymmTensorField kES_mean = (2.0/betaStar_) * (k_/(omega_ + low_omega_)) * S;
    volSymmTensorField sum_S = S;
    volScalarField k_mean = k_;
    volScalarField omega_mean = omega_;

    // Calcolo diagonali S^2 vettorializzato
    volScalarField squaredS_mean = sqr(S.component(symmTensor::XX)) 
                                 + sqr(S.component(symmTensor::YY)) 
                                 + sqr(S.component(symmTensor::ZZ));

    for (int i = 0; i < finestra_ - 1; ++i)
    {
        sum_uu += symm(UPrev_[i] * UPrev_[i]);
        sum_u += UPrev_[i];
        kES_mean += (2.0/betaStar_)*(kPrev_[i]/(omegaPrev_[i] + low_omega_)) * SPrev_[i];
        sum_S += SPrev_[i];
        k_mean += kPrev_[i];
        omega_mean += omegaPrev_[i];
        
        squaredS_mean += sqr(SPrev_[i].component(symmTensor::XX))
                       + sqr(SPrev_[i].component(symmTensor::YY))
                       + sqr(SPrev_[i].component(symmTensor::ZZ));
    }

    sum_u /= finestra_;
    volSymmTensorField L = (sum_uu/finestra_) - symm(sum_u * sum_u);
    kES_mean /= finestra_;
    volSymmTensorField S_avg = sum_S / finestra_;
    
    volScalarField k_mean_T = (k_mean/finestra_) + 0.5*(tr(sum_uu/finestra_) - (sum_u & sum_u)); 
    omega_mean /= finestra_;
    squaredS_mean /= finestra_;

    volScalarField SSquared_mean = sqr(S_avg.component(symmTensor::XX))
                                 + sqr(S_avg.component(symmTensor::YY))
                                 + sqr(S_avg.component(symmTensor::ZZ));

    omega_mean = (1.0 / (k_mean_T + low_k_)) * ((k_mean/finestra_ * omega_mean) + (this->nu() / betaStar_) * (squaredS_mean - SSquared_mean));
    volSymmTensorField Z = kES_mean - (2.0/betaStar_) * (k_mean_T / (omega_mean + low_omega_)) * S_avg;

    dynamicCmu_ = (Z && L) / ((Z && Z) + low_Z_);
    
    if (this->runTime_.timeIndex() < 10) dynamicCmu_ = 0.09;
    
    // Clipping vettorializzato
    dynamicCmu_ = min(max(dynamicCmu_, dynamicCmuMin_.value()), dynamicCmuMax_.value());

    // Aggiornamento Cronologia
    for (int i = 0; i < finestra_ - 2; ++i)
    {
        UPrev_[i] = UPrev_[i+1];
        kPrev_[i] = kPrev_[i+1];
        omegaPrev_[i] = omegaPrev_[i+1];
        SPrev_[i] = SPrev_[i+1];
    }

    UPrev_[finestra_ - 2] = this->U_;
    kPrev_[finestra_ - 2] = this->k_;
    omegaPrev_[finestra_ - 2] = this->omega_;
    SPrev_[finestra_ - 2] = S;

    correctNut();
}

} // End namespace RASModels
} // End namespace Foam
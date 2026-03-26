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

Class
    Foam::RASModels::kOmegaDyn

Description
    DYNAMIC k-omega turbulence model for incompressible flows.
    Default dynamic coefficient (dynamicCmu) clipping procedure -> [0;0.2]

    The default model coefficients are
    \verbatim
        kOmegaDynCoeffs
        {
            betaStar    0.09;
            gamma       0.52;
            beta        0.072;
            alphak      0.5;
            alphaOmega  0.5;
            dynamicCmuMin  0.0;
            dynamicCmuMax  0.2;
            nPrecedents    2;
        }
    \endverbatim

SourceFiles
    kOmegaDyn.C

\*---------------------------------------------------------------------------*/
#ifndef kOmegaDyn_H
#define kOmegaDyn_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                       Class kOmegaDyn Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class kOmegaDyn
:
    public eddyViscosity<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients
        dimensionedScalar betaStar_;
        dimensionedScalar beta_;
        dimensionedScalar gamma_;
        dimensionedScalar alphaK_;
        dimensionedScalar alphaOmega_;

        // dynamicCmu clipping limits
        dimensionedScalar dynamicCmuMin_;
        dimensionedScalar dynamicCmuMax_;

        const int finestra_;
        IOdictionary physicalProperties_;

        dimensionedScalar low_omega_;
        dimensionedScalar low_k_;
        dimensionedScalar low_Z_;

        // Fields
        volScalarField k_;
        volScalarField omega_;

        volSymmTensorField S;
        volTensorField gU;
        volScalarField dynamicCmu_;

        // Dynamic Procedure Buffer
        PtrList<volVectorField> UPrev_;
        PtrList<volScalarField> kPrev_;
        PtrList<volScalarField> omegaPrev_;
        PtrList<volSymmTensorField> SPrev_; // Added line
        PtrList<volVectorField> Uprint_;
        PtrList<volScalarField> kprint_;
        PtrList<volScalarField> omegaprint_;
        PtrList<volTensorField> gUprint_;

    // Private Member Functions
        virtual void correctNut();
        virtual tmp<fvScalarMatrix> kSource() const;
        virtual tmp<fvScalarMatrix> omegaSource() const;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    // Runtime type info
    TypeName("kOmegaDyn");

    // Constructors

        kOmegaDyn
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        // Do not allow copy
        kOmegaDyn(const kOmegaDyn&) = delete;

    // Destructor
    virtual ~kOmegaDyn() {}

    // Member Functions

        // Read coeff.
        virtual bool read();

        // Return effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return volScalarField::New("DkEff", alphaK_*this->nut_ + this->nu());
        }
        
        // Return effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return volScalarField::New("DomegaEff", alphaOmega_*this->nut_ + this->nu());
        }

        // Return turbulent kinetic energy
        virtual tmp<volScalarField> k() const { return k_; }

        // Return specific dissipation rate
        virtual tmp<volScalarField> omega() const { return omega_; }

        // Return turbulence dissipation
        virtual tmp<volScalarField> epsilon() const                     
        {
            return volScalarField::New
            (
                "epsilon",
                k_*omega_,
                omega_.boundaryField().types()
            );
        }

        // Solve equations and correct nut
        virtual void correct();

    // Member Operators
        // disallow assignment
        void operator=(const kOmegaDyn&) = delete;
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace RASModels
} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    # include "kOmegaDyn.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
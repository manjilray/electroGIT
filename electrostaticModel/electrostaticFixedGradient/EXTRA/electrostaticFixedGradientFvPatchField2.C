/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "electrostaticFixedGradientFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::electrostaticFixedGradientFvPatchField<Type>::electrostaticFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    gradient_(p.size(), Zero)
{}


template<class Type>
Foam::electrostaticFixedGradientFvPatchField<Type>::electrostaticFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    gradient_("gradient", dict, p.size())
{
    evaluate();
}


template<class Type>
Foam::electrostaticFixedGradientFvPatchField<Type>::electrostaticFixedGradientFvPatchField
(
    const electrostaticFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    gradient_(ptf.gradient_, mapper)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::electrostaticFixedGradientFvPatchField<Type>::electrostaticFixedGradientFvPatchField
(
    const electrostaticFixedGradientFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    gradient_(ptf.gradient_)
{}


template<class Type>
Foam::electrostaticFixedGradientFvPatchField<Type>::electrostaticFixedGradientFvPatchField
(
    const electrostaticFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::electrostaticFixedGradientFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    gradient_.autoMap(m);
}


template<class Type>
void Foam::electrostaticFixedGradientFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const electrostaticFixedGradientFvPatchField<Type>& fgptf =
        refCast<const electrostaticFixedGradientFvPatchField<Type>>(ptf);

    gradient_.rmap(fgptf.gradient_, addr);
}


template<class Type>
void Foam::electrostaticFixedGradientFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::electrostaticFixedGradientFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>(new Field<Type>(this->size(), pTraits<Type>::one));
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::electrostaticFixedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::electrostaticFixedGradientFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::electrostaticFixedGradientFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return gradient();
}


template<class Type>
void Foam::electrostaticFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    gradient_.writeEntry("gradient", os);
}


// ************************************************************************* //

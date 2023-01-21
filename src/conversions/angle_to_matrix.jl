# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related to the conversion from Euler angles to DCM.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export angle_to_matrix

# Mimics DCM constructor so code does not need to be touched beyond s/DCM/mkrotmat/g
function mkrotmat(m11,m12,m13,m21,m22,m23,m31,m32,m33)
    [m11 m12 m13; m21 m22 m23; m31 m32 m33]
end

"""
    angle_to_matrix(θ₁::Number, θ₂::Number, θ₃::Number, rot_seq::Symbol = :ZYX)
    angle_to_matrix(Θ::EulerAngles)

Convert the Euler angles `θ₁`, `θ₂`, and `θ₃` [rad] with the rotation sequence
`rot_seq` to a direction cosine matrix.

The input values of the origin Euler angles can also be passed inside the
structure `Θ` (see [`EulerAngles`](@ref)).

The rotation sequence is defined by a `:Symbol`. The possible values are:
`:XYX`, `XYZ`, `:XZX`, `:XZY`, `:YXY`, `:YXZ`, `:YZX`, `:YZY`, `:ZXY`, `:ZXZ`,
`:ZYX`, and `:ZYZ`. If no value is specified, then it defaults to `:ZYX`.

# Remarks

This function assigns `dcm = A3 * A2 * A1` in which `Ai` is the DCM related with
the *i*-th rotation, `i Є [1,2,3]`.

# Example

```jldoctest
julia> dcm = angle_to_matrix(pi / 2, pi / 3, pi / 4, :ZYX)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  3.06162e-17  0.5       -0.866025
 -0.707107     0.612372   0.353553
  0.707107     0.612372   0.353553
```
"""
function angle_to_matrix(θ₁::Number, θ₂::Number, θ₃::Number, rot_seq::Symbol = :ZYX)
    # Compute the sines and cosines.
    s₁, c₁ = sincos(θ₁)
    s₂, c₂ = sincos(θ₂)
    s₃, c₃ = sincos(θ₃)

    if rot_seq == :ZYX
        return mkrotmat(
                 c₂ * c₁,                c₂ * s₁,             -s₂ ,
            s₃ * s₂ * c₁ - c₃ * s₁, s₃ * s₂ * s₁ + c₃ * c₁, s₃ * c₂,
            c₃ * s₂ * c₁ + s₃ * s₁, c₃ * s₂ * s₁ - s₃ * c₁, c₃ * c₂
        )'
    elseif rot_seq == :XYX
        return mkrotmat(
              c₂,               s₁ * s₂,               -c₁ * s₂,
            s₂ * s₃, -s₁ * c₂ * s₃ + c₁ * c₃, c₁ * c₂ * s₃ + s₁ * c₃,
            s₂ * c₃, -s₁ * c₃ * c₂ - c₁ * s₃, c₁ * c₃ * c₂ - s₁ * s₃
        )'
    elseif rot_seq == :XYZ
        return mkrotmat(
             c₂ * c₃,  s₁ * s₂ * c₃ + c₁ * s₃, -c₁ * s₂ * c₃ + s₁ * s₃,
            -c₂ * s₃, -s₁ * s₂ * s₃ + c₁ * c₃,  c₁ * s₂ * s₃ + s₁ * c₃,
                s₂,             -s₁ * c₂,                 c₁ * c₂
        )'
    elseif rot_seq == :XZX
        return mkrotmat(
               c₂,               c₁ * s₂,                 s₁ * s₂,
            -s₂ * c₃,  c₁ * c₃ * c₂ - s₁ * s₃,  s₁ * c₃ * c₂ + c₁ * s₃,
             s₂ * s₃, -c₁ * c₂ * s₃ - s₁ * c₃, -s₁ * c₂ * s₃ + c₁ * c₃
        )'

    elseif rot_seq == :XZY
        return mkrotmat(
            c₃ * c₂, c₁ * c₃ * s₂ + s₁ * s₃, s₁ * c₃ * s₂ - c₁ * s₃,
              -s₂,             c₁ * c₂,                s₁ * c₂,
            s₃ * c₂, c₁ * s₂ * s₃ - s₁ * c₃, s₁ * s₂ * s₃ + c₁ * c₃
        )'
    elseif rot_seq == :YXY
        return mkrotmat(
            -s₁ * c₂ * s₃ + c₁ * c₃,  s₂ * s₃, -c₁ * c₂ * s₃ - s₁ * c₃,
                       s₁ * s₂,         c₂,              c₁ * s₂,
             s₁ * c₃ * c₂ + c₁ * s₃, -s₂ * c₃,  c₁ * c₃ * c₂ - s₁ * s₃
        )'
    elseif rot_seq == :YXZ
        return mkrotmat(
             c₁ * c₃ + s₂ * s₁ * s₃, c₂ * s₃, -s₁ * c₃ + s₂ * c₁ * s₃,
            -c₁ * s₃ + s₂ * s₁ * c₃, c₂ * c₃,  s₁ * s₃ + s₂ * c₁ * c₃,
                  s₁ * c₂,             -s₂,         c₂ * c₁
        )'
    elseif rot_seq == :YZX
        return mkrotmat(
                       c₁ * c₂,         s₂,              -s₁ * c₂,
            -c₃ * c₁ * s₂ + s₃ * s₁,  c₂ * c₃,  c₃ * s₁ * s₂ + s₃ * c₁,
             s₃ * c₁ * s₂ + c₃ * s₁, -s₃ * c₂, -s₃ * s₁ * s₂ + c₃ * c₁
        )'
    elseif rot_seq == :YZY
        return mkrotmat(
            c₁ * c₃ * c₂ - s₁ * s₃, s₂ * c₃, -s₁ * c₃ * c₂ - c₁ * s₃,
                -c₁ * s₂,             c₂,               s₁ * s₂,
            c₁ * c₂ * s₃ + s₁ * c₃, s₂ * s₃, -s₁ * c₂ * s₃ + c₁ * c₃
        )'
    elseif rot_seq == :ZXY
        return mkrotmat(
            c₃ * c₁ - s₂ * s₃ * s₁, c₃ * s₁ + s₂ * s₃ * c₁, -s₃ * c₂,
                -c₂ * s₁,                c₂ * c₁,              s₂,
            s₃ * c₁ + s₂ * c₃ * s₁, s₃ * s₁ - s₂ * c₃ * c₁,  c₂ * c₃
        )'
    elseif rot_seq == :ZXZ
        return mkrotmat(
            -s₁ * c₂ * s₃ + c₁ * c₃, c₁ * c₂ * s₃ + s₁ * c₃, s₂ * s₃,
            -s₁ * c₃ * c₂ - c₁ * s₃, c₁ * c₃ * c₂ - s₁ * s₃, s₂ * c₃,
                       s₁ * s₂,               -c₁ * s₂,         c₂
        )'
    elseif rot_seq == :ZYZ
        return mkrotmat(
             c₁ * c₃ * c₂ - s₁ * s₃,  s₁ * c₃ * c₂ + c₁ * s₃, -s₂ * c₃,
            -c₁ * c₂ * s₃ - s₁ * c₃, -s₁ * c₂ * s₃ + c₁ * c₃,  s₂ * s₃,
                       c₁ * s₂,                 s₁ * s₂,          c₂
        )'
    else
        throw(ArgumentError("The rotation sequence :$rot_seq is not valid."))
    end
end

angle_to_matrix(Θ::EulerAngles) = angle_to_matrix(Θ.a1, Θ.a2, Θ.a3, Θ.rot_seq)

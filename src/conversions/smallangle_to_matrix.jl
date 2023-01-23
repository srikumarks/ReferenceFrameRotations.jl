# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related to the conversion from small Euler angles to Matrix.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export smallangle_to_matrix

"""
    smallangle_to_matrix(θx::Number, θy::Number, θz::Number; normalize = true)

Create a direction cosine matrix from three small rotations of angles `θx`,
`θy`, and `θz` [rad] about the axes X, Y, and Z, respectively.

If the keyword `normalize` is `true`, then the matrix will be normalized using
the function `orthonormalize`.

# Example

```jldoctest
julia> smallangle_to_matrix(+0.01, -0.01, +0.01)
Matrix{Float64}:
  0.9999     0.00989903  0.010098
 -0.009999   0.999901    0.00989802
 -0.009999  -0.009998    0.9999

julia> smallangle_to_matrix(+0.01, -0.01, +0.01; normalize = false)
Matrix{Float64}:
  1.0    0.01  0.01
 -0.01   1.0   0.01
 -0.01  -0.01  1.0
```
"""
@inline function smallangle_to_matrix(
    θx::T1,
    θy::T2,
    θz::T3;
    normalize = true
) where {T1<:Number, T2<:Number, T3<:Number}
    # Since we might orthonormalize `D`, we need to get the float to avoid type
    # instabilities.
    T = float(promote_type(T1, T2, T3))

    D = mkrotmat(
          1, +θz, -θy,
        -θz,   1, +θx,
        +θy, -θx,   1
    )

    if normalize
        return orthonormalize(D)
    else
        return D
    end
end

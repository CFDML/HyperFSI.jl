using TestItemRunner

   """
    are_cols_approx_equal(A, B; atol=1e-8)

    Check if matrices `A` and `B` contain the same columns (order-independent), 
    using floating-point tolerance `atol`.
    """
    function are_cols_approx_equal(A::AbstractMatrix, B::AbstractMatrix; atol::Real=1e-8)
        size(A) == size(B) || return false
        for col in eachcol(A)
            any(≈(col; atol), eachcol(B)) || return false
        end
        return true
    end

# Only run tests that are not tagged with :skipci due to performance issues in CI
@run_package_tests verbose=true filter=ti->!(:skipci in ti.tags)
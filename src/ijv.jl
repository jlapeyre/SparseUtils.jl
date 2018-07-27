module IJV

function dropzeros!(I, J, V::AbstractArray{T}) where T
    mask = V .== zero(T)
    map(x -> deleteat!(x, mask), (I, J, V))
end
dropzeros(I, J, V::AbstractArray) = dropzeros!(copy(I), copy(J), copy(V))

end

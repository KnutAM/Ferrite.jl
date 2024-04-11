vertex_coordinates(RF::Type{<:AbstractRefShape}) = vertex_coordinates(Float64, RF)
vertex_coordinates(::Type{T}, ::Type{RefLine}) where {T<:Number} = (Vec{1,T}((-1,)), Vec{1,T}((1,)))
function vertex_coordinates(::Type{T}, ::Type{RefTriangle}) where {T<:Number}
    vf(args...) = Vec{2, T}(args)
    return (vf(1, 0), vf(0, 1), vf(0, 0))
end
function vertex_coordinates(::Type{T}, ::Type{RefQuadrilateral}) where {T<:Number}
    vf(args...) = Vec{2, T}(args)
    return (vf(-1, -1), vf(1, -1), vf(1, 1), vf(-1, 1))
end
function vertex_coordinates(::Type{T}, ::Type{RefTetrahedron}) where {T<:Number}
    vf(args...) = Vec{3, T}(args)
    return (vf(0, 0, 0), vf(1, 0, 0), vf(0, 1, 0), vf(0, 0, 1))
end
function vertex_coordinates(::Type{T}, ::Type{RefHexahedron}) where {T<:Number}
    vf(args...) = Vec{3, T}(args)
    return (vf(-1, -1, -1), vf(1, -1, -1), vf(1, 1, -1), vf(-1, 1, -1), vf(-1, -1, 1), vf(1, -1, 1), vf(1, 1, 1), vf(-1, 1, 1))
end
function vertex_coordinates(::Type{T}, ::Type{RefPrism}) where {T<:Number}
    vf(args...) = Vec{3, T}(args)
    return (vf(0, 0, 0), vf(1, 0, 0), vf(0, 1, 0), vf(0, 0, 1), vf(1, 0, 1), vf(0, 1, 1))
end

#=
vctest(::Type{RF}) where RF = reference_coordinates(Lagrange{RF,1}())
for RS in (RefLine, RefTriangle, RefQuadrilateral, RefHexahedron, RefPrism)
    println(RS, ": ", vctest(RS) ≈ collect(vertex_coordinates(RS)))
end
=#

refshape_edges(::Type{<:RefHexahedron}) = ((1,2), (2,3), (3,4), (4,1), (5,6), (6,7), (7,8), (8,5), (1,5), (2,6), (3,7), (4,8))
refshape_edges(::Type{<:RefPrism}) = ((2,1), (1,3), (1,4), (3,2), (2,5), (3,6), (4,5), (4,6), (6,5))

refshape_faces(::Type{<:RefQuadrilateral}) = ((1,2), (2,3), (3,4), (4,1))
refshape_faces(::Type{<:RefHexahedron}) = ((1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8))
refshape_faces(::Type{<:RefPrism}) = (1,3,2), (1,2,5,4), (3,1,4,6), (2,3,6,5), (4,5,6)

refshape_faces_by_edge(::Type{<:RefQuadrilateral}) = ntuple(_->(), 4)
# ordering not well-defined...
refshape_faces_by_edge(::Type{<:RefHexahedron}) = ((1,2,3,4), (1,9,5,10), (2,11,6,10), (3,12,7,11), (4,12,8,9), (5,6,7,8))
refshape_faces_by_edge(::Type{<:RefPrism}) = ((1,2,4), (1,3,7,5), (2,6,8,3), (4,6,9,5), (7,8,9))

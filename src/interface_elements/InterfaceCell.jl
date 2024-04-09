refshape_edges(::RefHexahedron) = ((1,2), (2,3), (3,4), (4,1), (5,6), (6,7), (7,8), (8,5), (1,5), (2,6), (3,7), (4,8))
refshape_edges(::RefPrism) = ((2,1), (1,3), (1,4), (3,2), (2,5), (3,6), (4,5), (4,6), (6,5))

refshape_faces(::RefQuadrilateral) = ((1,2), (2,3), (3,4), (4,1))
refshape_faces(::RefHexahedron) = ((1,4,3,2), (1,2,6,5), (2,3,7,6), (3,4,8,7), (1,5,8,4), (5,6,7,8))
refshape_faces(::RefPrism) = (1,3,2), (1,2,5,4), (3,1,4,6), (2,3,6,5), (4,5,6)

struct InterfaceInterpolation{RefShape, IPhere, IPthere} <: Interpolation{RefShape, Nothing, Nothing}
    here::IPhere
    there::IPthere
    function InterfaceInterpolation(ip_here::ScalarInterpolation{EmbeddedRefShape}, ip_there::ScalarInterpolation{EmbeddedRefShape}) where EmbeddedRefShape
        RefShape = interface_refshape_from_embedded(RefShape)
        return new{RefShape, typeof(ip_here), typeof(ip_there)}(ip_here, ip_there)
    end
    function InterfaceInterpolation(ip_here::VectorInterpolation{vdim, EmbeddedRefShape}, ip_there::ScalarInterpolation{vdim, EmbeddedRefShape}) where {vdim, EmbeddedRefShape}
        RefShape = interface_refshape_from_embedded(RefShape)
        return new{RefShape, typeof(ip_here), typeof(ip_there)}(ip_here, ip_there)
    end
end   

#interface_refshape_from_embedded(::Type{RefPoint}) = RefLine # not implemented
interface_refshape_from_embedded(::Type{RefLine}) = RefQuadrilateral
interface_refshape_from_embedded(::Type{RefTriangle}) = RefPrism
interface_refshape_from_embedded(::Type{RefQuadrilateral}) = RefHexahedron

function vertexdof_indices(ip::InterfaceInterpolation{RS}) where {RS <: Union{RefHexahedron, RefPrism}}
    dofs_here = vertexdof_indices(ip.here)
    nhere = sum(length, dofs_here)
    dofs_there = map(dofs -> map(d -> d + nhere, dofs), vertexdof_indices(ip.there))
    return tuple(dofs_here..., dofs_there...)
end
function vertexdof_indices(ip::InterfaceInterpolation{RefQuadrilateral})
    dofs_here = vertexdof_indices(ip.here)
    nhere = sum(length, dofs_here)
    dofs_there = map(dofs -> map(d -> d + nhere, dofs), reverse(vertexdof_indices(ip.there)))
    return tuple(dofs_here..., dofs_there...)
end

function edgedof_interior_indices_hereandthere(ip::InterfaceInterpolation{RS}) where {RS <: Union{RefHexahedron, RefPrism}}
    nv_here = sum(length, vertexdof_indices(ip.here))
    nv_there = sum(length, vertexdof_indices(ip.there))
    fdofs_i_here = facedof_interior_indices(ip.here)
    Δn_there = nv_here + sum(length, fdofs_i_here)
    dofs_here = map(dofs -> map(d -> d + nv_there, dofs), fdofs_i_here)
    dofs_there = map(dofs -> map(d -> d + Δn_there, dofs), facedof_interior_indices(ip.there))
    return dofs_here, dofs_there
end

function edgedof_interior_indices(ip::InterfaceInterpolation{RefHexahedron})
    dofs_here, dofs_there = edgedof_interior_indices_hereandthere(ip)
    return tuple(dofs_here..., dofs_there..., ntuple(_->(), 4)...)
end
function edgedof_interior_indices(ip::InterfaceInterpolation{RefPrism})
    dofs_here, dofs_there = edgedof_interior_indices_hereandthere(ip)
    return tuple(dofs_here[1:2]..., (), dofs_here[3], (), (), dofs_there...)
end

function edgedof_indices(ip::InterfaceInterpolation{RS}) where {RS <: AbstractRefShape{3}}
    dofs_int = edgedof_interior_indices(ip)
    dofs_vert = vertexdof_indices(ip)
    return map(dofs_int, refshape_edges(RS)) do (dints, vnrs)
        dv1 = dofs_vert[vnrs[1]] # dofs at vertex 1
        dv2 = dofs_vert[vnrs[2]] # dofs at vertex 2
        tuple(dv1..., dv2..., dints...)
    end
end

function facedof_interior_indices_hereandthere(ip::InterfaceInterpolation{RS}) where {RS <: AbstractRefShape{3}}
    dofs_here = celldof_interior_indices(ip.here)
    nve_here = getnbasefunctions(ip.here) - length(dofs_here)
    dofs_there = celldof_interior_indices(ip.there)
    nve_there = getnbasefunctions(ip.here) - length(dofs_here)
    Δn_here = nve_here + length(dofs_here)
    return (map(d -> d + nve_there, dofs_here),
            map(d -> d + Δn_here,   dofs_there))
end

function facedof_interior_indices(ip::InterfaceInterpolation{RefShape}) where RefShape
    dofs_here, dofs_there = facedof_interior_indices_hereandthere(ip)
    
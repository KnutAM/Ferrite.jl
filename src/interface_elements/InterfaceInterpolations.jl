combinetuples(args...) = reduce((a,b) -> tuple(a..., b...), args; init=())

# Support distribution of vertexdof_indices directly
#=
function vertexdof_indices(ip::VectorizedInterpolation{vdim}) where vdim
    return map(vertexdof_indices(ip.ip)) do d
=#        

struct InterfaceInterpolation{RefShape, IPhere, IPthere} <: Interpolation{RefShape, Nothing, Nothing}
    here::IPhere
    there::IPthere
    function InterfaceInterpolation(ip_here::ScalarInterpolation{EmbeddedRefShape}, ip_there::ScalarInterpolation{EmbeddedRefShape}) where EmbeddedRefShape
        RefShape = interface_refshape_from_embedded(EmbeddedRefShape)
        return new{RefShape, typeof(ip_here), typeof(ip_there)}(ip_here, ip_there)
    end
    function InterfaceInterpolation(ip_here::VectorInterpolation{vdim, EmbeddedRefShape}, ip_there::VectorInterpolation{vdim, EmbeddedRefShape}) where {vdim, EmbeddedRefShape}
        RefShape = interface_refshape_from_embedded(EmbeddedRefShape)
        return new{RefShape, typeof(ip_here), typeof(ip_there)}(ip_here, ip_there)
    end
end
n_components(ip::InterfaceInterpolation) = n_components(ip.here)
adjust_dofs_during_distribution(ip::InterfaceInterpolation) = adjust_dofs_during_distribution(ip.here)

#interface_refshape_from_embedded(::Type{RefPoint}) = RefLine # RefPoint not implemented
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
    dofs_there = map(dofs -> map(d -> d + nhere, dofs), #=reverse(=#vertexdof_indices(ip.there)#=)=#)
    return tuple(dofs_here..., dofs_there...)
end

function edgedof_interior_indices_hereandthere(ip::InterfaceInterpolation{RS}) where {RS <: Union{RefHexahedron, RefPrism}}
    nv_there = sum(length, vertexdof_indices(ip.there))
    nve_there = getnbasefunctions(ip.here) - length(celldof_interior_indices(ip.here))
    dofs_here  = map(dofs -> map(d -> d + nv_there,  dofs), facedof_interior_indices(ip.here))
    dofs_there = map(dofs -> map(d -> d + nve_there, dofs), facedof_interior_indices(ip.there))
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
        vdofs = combinetuples(i -> dofs_vert[i], vnrs)
        tuple(vdofs..., dints...)
    end
end

function facedof_interior_indices_hereandthere(ip::InterfaceInterpolation{RS}) where {RS <: AbstractRefShape{3}}
    dofs_here = celldof_interior_indices(ip.here)
    dofs_there = celldof_interior_indices(ip.there)
    nve_there = getnbasefunctions(ip.here) - length(dofs_here)
    n_here = getnbasefunctions(ip.here)
    # Cannot have more than 1 interior facedof in 3d
    @assert length(dofs_here) ≤ 1 ≥ length(dofs_there)
    return (map(d -> d + nve_there, dofs_here),
            map(d -> d + n_here,    dofs_there))
end

function facedof_interior_indices(ip::InterfaceInterpolation{RefHexahedron})
    dofs_here, dofs_there = facedof_interior_indices_hereandthere(ip)
    return tuple(dofs_here, ntuple(_->(), 4)..., dofs_there)
end

function facedof_interior_indices(ip::InterfaceInterpolation{RefPrism})
    dofs_here, dofs_there = facedof_interior_indices_hereandthere(ip)
    return tuple(dofs_here, ntuple(_->(), 3)..., dofs_there)
end

function facedof_interior_indices(ip::InterfaceInterpolation{RefQuadrilateral})
    dofs_here = celldof_interior_indices(ip.here)
    dofs_there = reverse(celldof_interior_indices(ip.there))
    n_here = getnbasefunctions(ip.here)
    nv_there = sum(length, vertexdof_indices(ip.there))
    return tuple(map(d -> d + nv_there, dofs_here), (), 
                 map(d -> d + n_here, dofs_there), ())
end

function facedof_indices(ip::InterfaceInterpolation{RefShape}) where RefShape
    dofs_int = facedof_interior_indices(ip)
    dofs_edge = edgedof_interior_indices(ip)
    dofs_vert = vertexdof_indices(ip)
    return map(dofs_int, refshape_faces(RefShape), refshape_faces_by_edge(RefShape)) do (dints, vnrs, enrs)
        vdofs = combinetuples(map(i -> dofs_vert[i], vnrs))
        edofs = combinetuples(map(i -> dofs_edge[i], enrs))
        tuple(vdofs..., edofs..., dints...)
    end
end

using Ferrite: nnodes_per_cell

abstract type Load{dim} end
struct NodalLoad{dim} <: Load{dim}
    nodeset_name::AbstractString
    F::NTuple{dim}
end
struct LinearLoad{dim} <: Load{dim}
    faceset_name::AbstractString
    F::NTuple{dim}
end

struct FEModel{dim}
    mat::Material{dim}
    grid::Grid
    elemvol::Vector{<:Real}

    constraints::Vector{Dirichlet}
    loads::Vector{<:Load{dim}}

    cellvalues::CellValues
    facevalues::FaceValues
    dh::DofHandler
    ch::ConstraintHandler
end

function FEModel(; mat::Material{dim}, grid::Grid, ip::Interpolation, qr::QuadratureRule{dim,shape}, constraints::Vector{Dirichlet}, loads::Vector{<:Load}) where {dim,shape}
    # element type and quadrature rule
    cellvalues = CellVectorValues(qr, ip)

    qr_order = length(getpoints(qr))
    face_qr = QuadratureRule{dim - 1,shape}(qr_order)
    facevalues = FaceVectorValues(face_qr, ip)

    # degrees of freedom
    dh = DofHandler(grid)
    add!(dh, :u, dim)
    close!(dh)

    # constraints
    ch = ConstraintHandler(dh)
    for cc in constraints
        add!(ch, cc)
    end
    close!(ch)

    # elemental volume vector
    elemvol = zeros(getncells(grid))
    for cell in CellIterator(dh)
        reinit!(cellvalues, cell)
        for q_point in 1:getnquadpoints(cellvalues)
            elemvol[cellid(cell)] += getdetJdV(cellvalues, q_point)
        end
    end

    return FEModel(mat, grid, elemvol, constraints, loads, cellvalues, facevalues, dh, ch)
end

function get_centers(model::FEModel{dim}) where {dim}
    centers = zeros(getncells(model.grid), dim)
    for cell in CellIterator(model.dh)
        id = cellid(cell)
        for node in getcoordinates(model.grid, id)
            centers[id, :] .+= node
        end
        centers[id, :] ./= nnodes_per_cell(model.grid, id)
    end
    return centers
end

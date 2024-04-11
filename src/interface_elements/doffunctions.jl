
"Dofs associated with 0-dimensional entities (Ciarlet: point)"
vertexdof_indices # Same as currently

"Dofs associated with 1-dimensional entities (Ciarlet: edge)"
curvedof_indices

"Dofs associated with 2-dimensional entities (Ciarlet: face)"
surfacedof_indices

"Dofs associated with 3-dimensional entities (Ciarlet: volume)"
volumedof_indices

"Dofs associated with codimension 2 entities (only applicable to 3d and 2d objects) (Ciarlet: ridges)"
edgedof_indices # ref meaning of edges(::AbstractCell)

"Dofs associated with codimension 1 entities (Ciarlet: facets)"
facedof_indices # ref meaning of faces(::AbstractCell)

"Dofs associated with codimension 0 entities (Ciarlet: cell)"
celldof_indices

# dof-distribution should go `vertex - curve - surface - volume`

# edges
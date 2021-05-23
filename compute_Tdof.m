function [Td] = compute_Tdof(Beam_elements, n_dof, n_nod_b, Tbeams)

Td = zeros(Beam_elements,n_nod_b*n_dof);
for e=1:Beam_elements
    for a=1:n_nod_b
        Td(e,a*n_dof - (n_dof-1):a*n_dof) = Tbeams(e,a)*n_dof-(n_dof-1):Tbeams(e,a)*n_dof;
    end
end

end
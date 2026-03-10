function mat_vec = get_material_property(mat_type)
% -------------------------------------------------------------------------
% Helper routine
%   get_material_property(mat_type)
%
% Returns a 1×3 row vector [ρ  E  ν] corresponding to the material name
% passed in 'mat_type'.  Add new cases here whenever you introduce a new
% material to your lamination schedule.
% -------------------------------------------------------------------------

    switch mat_type

        case 'Composite_1'
            E11 = 40e9;
            E22 = 1e9;
            G12 = 0.6e9;
            G23 = 0.5e9;
            G31 = G12;
            nu12 = 0.25;
            nu21 = E22*nu12/E11;
            rho = 1000;

        % --------------- Default / unknown material ----------------------
        otherwise
            error('Material "%s" not recognised in get_material_property.', mat_type);
    end

    % Return the triplet as a row-vector
    mat_vec = [rho, E11, E22, G12, G31, G23, nu12, nu21];
end
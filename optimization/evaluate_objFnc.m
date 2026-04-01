function obj = evaluate_objFnc(baseline_plate,baseline_stiff,dcp)
%
[Nurbs2D_plate,Nurbs2D_stiff] = update_baseline(baseline_plate,baseline_stiff,dcp);
%
semOpt = semOptMesh(Nurbs2D,N,shell_dof);
%

end
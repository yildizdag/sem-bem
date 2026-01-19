function k_loc = localStiffness1D(sem1D,el)
J = (sem1D.Q1)*sem1D.nodes(sem1D.conn(el,:),1);
k_loc = (sem1D.EA).*((sem1D.Q1./J)'*sem1D.V*(sem1D.Q1./J)).*J;
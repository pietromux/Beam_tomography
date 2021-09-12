function [phsp_f, R] = matrix_quad_transport(phsp, quadvalues)

% This function takes as inputs the initial 4D distribution and 3 quad current values and runs them thru matrix transport until the screen after the box
gammaBeta = mean(phsp(:,6));
dbs = 1.1;    % to be checked with Youna - should be consistent with GPT screen location

[ R , R_all, z_all ,pos] = transport_matrix_ff_focusing_green_var(quadvalues,gammaBeta,dbs);
phsp_f = phsp;
phsp_f(:,1:4) = (R*phsp(:,1:4)')';
end

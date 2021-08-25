function phsp_out = momentum_reassign(phsp, idx, r0p, flagdist)

phsp_out = phsp;
nspawn = sum(idx);
nfix = sum(~idx);
npart = nfix+nspawn;


%Set random angular distribution
if (strcmp(flagdist,'gaussian'))    
    phsp_out(idx,2) = r0p*randn(nspawn,1);
    phsp_out(idx,4) = r0p*randn(nspawn,1);
elseif (strcmp(flagdist,'uniform'))   
    r2p = rand(nspawn,1)*r0p*r0p;
    phip = rand(nspawn,1)*2*pi;
    phsp_out(idx,2)=sqrt(r2p).*cos(phip);
    phsp_out(idx,4)=sqrt(r2p).*sin(phip);
elseif (strcmp(flagdist,'fit'))
    s_attempt = 4*cov(phsp(~idx(:),1:4));
    phsp_attempt = mvnrnd([0 0 0 0],s_attempt,nspawn);
%    phsp_out(idx,2) = phsp_attempt(:,2);
%    phsp_out(idx,4) = phsp_attempt(:,4);
    phsp_out(idx,1:4) = phsp_attempt(:,1:4);
elseif (strcmp(flagdist,'randomwalk'))
    phsp_out(idx,2) = phsp_out(idx,2)+r0p*rand(1)*cos(rand(1)*2*pi);
    phsp_out(idx,4) = phsp_out(idx,4)+r0p*rand(1)*sin(rand(1)*2*pi);
end

end

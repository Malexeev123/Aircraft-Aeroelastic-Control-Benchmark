function idx = buildIndexStruct(Nm,Na2)
    % Utility: called once in SimRunner ctor to pre‑compute ranges.
    idx.q1  = 1:Nm;
    idx.q2  = Nm + (1:Nm);
    idx.qxi = 2*Nm + (1:Nm+1);
    idx.qGam= (3*Nm+1) + (1:Na2);
    idx.chi = (3*Nm+1+Na2) + (1:3);
end
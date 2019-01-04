function [F,G,H] = OraclePH(qc,ind)
    q = q0 + B * qc;
    if ind == 2 | ind == 4 | ind == 7 then
        F = (1/3) * q' * (r .* q .* abs(q)) + pr' * (Ar * q);
    end
    if ind == 3 | ind == 4 | ind == 6 | ind == 7 then
        G = B' * (r .* q .* abs(q) + Ar' * pr)
    end
    if ind == 5 | ind == 6 | ind == 7 then
        H = 2 * B' * diag([r .* abs(q)]) * B
    end
endfunction

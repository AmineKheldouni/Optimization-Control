
function [F,G] = OraclePG(qc,ind)
  q = q0+B*qc;
  if (ind == 2 | ind == 4) then
      F = (1/3) * (q)' * (r .* (q0 + B * qc) .* abs(q)) + pr' * (Ar*(q));
  end
  if (ind == 3 | ind == 4) then
      G = B' * (r .* (q) .* abs((q)) + Ar' * pr);
  end
endfunction

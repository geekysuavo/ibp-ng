
% load the input data file.
dat = load(argv(){1});
n = rows(dat);

% get the phi and psi matrices.
Phi = dat(:, 2 : 2 : end);
Psi = dat(:, 3 : 2 : end);

% get the bounds on phi.
phi_min = min(Phi, [], 2);
phi_max = max(Phi, [], 2);

% get the bounds on psi.
psi_min = min(Psi, [], 2);
psi_max = max(Psi, [], 2);

% loop over the residues.
for resid = 1 : n
  % comment.
  printf('{* resid %d *}\n', resid);

  % get residue ids.
  I = resid - 1;
  J = resid;
  K = resid + 1;

  % phi.
  if (resid > 1)
    phi = (phi_max(resid) + phi_min(resid)) / 2;
    dphi = phi_max(resid) - phi_min(resid);

    r_dihe(I, 'C', J, 'N', J, 'CA', J, 'C', phi, dphi);
  end

  % psi, omega.
  if (resid < n)
    psi = (psi_max(resid) + psi_min(resid)) / 2;
    dpsi = psi_max(resid) - psi_min(resid);

    r_dihe(J,  'N', J, 'CA', J, 'C', K,  'N', psi, dpsi);
    r_dihe(J, 'CA', J,  'C', K, 'N', K, 'CA', 180.0, 0.0);
  end
  printf('\n');
end

% add base restraints to keep ibp from branching on HA and
% terminal O, O2.
printf('{* base restraints *}\n');
for resid = 3 : n
  r_dihe(resid, 'N', resid, 'CA', resid, 'C', resid, 'HA', -119, 0);
end
r_dihe(n, 'HA', n, 'C', n, 'CA', n, 'O', 128.5, 0);
r_dihe(n, 'CA', n, 'O', n, 'C', n, 'O2', 171.5, 0);


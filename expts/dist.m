
% set the expansion value.
expand = 1.0;

% get the input data files.
files = glob(argv(){1});
nfiles = length(files);

% build the distance matrices.
D = cell(nfiles, 1);
for ii = 1 : nfiles
  % load the input file.
  dat = load(files{ii});

  % get the x, y, and z vectors.
  id = dat(:, 1);
  x = dat(:, 2);
  y = dat(:, 3);
  z = dat(:, 4);

  % compute the distance matrix.
  D{ii} = sqrt(bsxfun(@minus, x, x').^2 + ...
               bsxfun(@minus, y, y').^2 + ...
               bsxfun(@minus, z, z').^2);
end

% build an array of the distance matrices.
dmat = zeros(rows(D{1}), columns(D{1}), length(D));
for ii = 1 : length(D)
  dmat(:, :, ii) = D{ii};
end
D = dmat;

% compute the distance lower and upper bounds.
Dmin = min(D, [], 3);
Dmax = max(D, [], 3);
Dmean = mean(D, 3);

% find all proton pairs having a distance within 6A
% in one or more of the structures.
[i, j, d] = find(Dmin <= 6);

% loop over the identified pairs.
printf('\n{* distance restraints *}\n');
for ii = 1 : length(d)
  if (i(ii) < j(ii) && abs(i(ii) - j(ii)) > 5)
    % extract the bounds and the mean.
    L = Dmin(i(ii), j(ii));
    U = Dmax(i(ii), j(ii));
    M = Dmean(i(ii), j(ii));

    % output the restraint.
    r_dist(i(ii), 'H1', j(ii), 'H1', M, M - L + expand, U - M + expand);
  end
end


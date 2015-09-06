function y = sosfiltfilt(sos, x)
% SOSFILTFILT Zero-phase filtering with second-order sections
%
% SOSFILTFILT implements zero-phase forward and reverse digital filtering
% using the second-order sections filter implementation.  The resulting
% filter has zero phase distortion, but the squared magnitude of the original
% filter.
%
% usage: y = sosfiltfilt(sos, x);
%
%   sos   second-order sections model for filter
%   x     signal to be filter
%
%   y     resulting filtered signal
%
% See also SOSFILT and FILTFILT.

% Shourov K. Chatterji
% shourov@ligo.mit.edu
% 2003-Jul-14

% parse command line arguments
error(nargchk(2,2,nargin));

% force row vector
if size(x, 1) > size(x, 2),
  iscolumn = 1;
else
  iscolumn = 0;
end
x = x(:)';

% get zero-pole-gain representation
[z, p, k] = sos2zp(sos);

% approximate duration of startup transient
transient = ceil(max(pi ./ (abs(p) .* angle(p))));

% limit transient to data length
transient = min(transient, length(x) - 1);

% pad with reflected data to reduce startup transients
prepend = 2 * x(1) - x(transient+1:-1:2);
append = 2 * x(end) - x(end-1:-1:end-transient);
x = [prepend x append];

% forward filter data
y = sosfilt(sos, x);

% reverse result
y = fliplr(y);

% reverse filter data
y = sosfilt(sos, y);

% reverse result
y = fliplr(y);

% remove transients
y = y(transient+1:end-transient);

% return same vector orientation
if iscolumn == 1,
  y = y(:);
end
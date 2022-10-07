function E = wgs84Ellipsoid(lengthUnit)
%% wgs84Ellipsoid  generate a WGS84 referenceEllipsoid
%
%%%  Inputs
%
% * lengthUnit: currently not implemented, for compatibility either nargin=0 or 'm' or 'meter' is accepted without conversion
%
%%% Outputs
%
% * E: referenceEllipsoid
arguments
  lengthUnit (1,1) string = "m"
end
%% get ellipsoid
E = matmap3d.referenceEllipsoid('wgs84', lengthUnit);

end
%%
%  Copyright (c) 2014-2018 Michael Hirsch, Ph.D.
%
%  Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
%  1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
%  2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

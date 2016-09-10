% SphericalHarmonics.m %
% http://mathworld.wolfram.com/SphericalHarmonic.html %

function value = SphericalHarmonics(L, M, theta, phi)
  YLMTemp = legendre(L, cos(theta), 'norm');
  factor = (-1)^abs(M);
  value = factor .* YLMTemp(abs(M)+1, :) .* exp(j * abs(M) .* phi') / sqrt(2 * pi);
  if M < 0
    value = factor * conj(value);
  elseif L==0
    value = value / sqrt(2); 
  end;
  
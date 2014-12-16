function fx = cauchyPDF2(x, x0, gamma)

fx = ( 1./ (pi*gamma .* (1 + ((x-x0)./gamma).^2 ) ) );
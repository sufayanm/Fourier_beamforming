function kmig = find_kmig(kz, kx, ku)
kv2= (kx-ku).^2 ;
ku2 = ku.^2 ;
kz2 = kz.^2 ;

% kmig = sqrt(kz2.^2 + 2*(ku2 + kv2).*kz2 + ku2.^2 + kv2.^2 - 2*ku2.*kv2)./(2*kz) ;
kmig = kz2/4 + 0.5*(ku2 + kv2 ) + ((ku2- kv2).^2)./(4.*kz2) ;

kmig(isinf(kmig)|isnan(kmig)) = 0 ;
kmig = sqrt(kmig) ;

kmig(kz<0 , : ,:) = -kmig(kz<0 , : ,:) ;
% kmig(kz<0 , : ,:) = 0 ;
end
function stdev_ro = sum_ro(trenutni_pravac_brzine, trenutno_teta, donja_granica, gornja_granica) % Generisanje šuma za poteg

trenutno_teta = trenutno_teta*180/pi;
trenutno_teta = mod(trenutno_teta,360);

trenutni_pravac_brzine = trenutni_pravac_brzine*180/pi;
trenutni_pravac_brzine = mod(trenutni_pravac_brzine,360);

if (abs(trenutno_teta-trenutni_pravac_brzine) <= 90)
    rad = abs(trenutni_pravac_brzine - trenutno_teta);
elseif (abs(trenutno_teta-trenutni_pravac_brzine) <= 180)
    rad = 180 - abs(trenutni_pravac_brzine - trenutno_teta);
elseif (abs(trenutno_teta-trenutni_pravac_brzine) <= 270)
    rad = abs(trenutni_pravac_brzine - trenutno_teta) - 180;
else
    rad = 360 - abs(trenutni_pravac_brzine - trenutno_teta);
end

stdev_ro = -donja_granica*rad/90 + gornja_granica;


    

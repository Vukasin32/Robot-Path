clear
close all
clc

%% Inicijalizacija

podaci = csvread('observations.csv'); % Preskačemo prvi red (zaglavlje)
ro_izmereno = podaci(:,1);
teta_izmereno = podaci(:,2);

broj_prikazanih_cestica = 100;

N = 2000;
cestice = zeros(N,3); % 3 komp., x i y koordinata i ugao brzine 

for i = 1:N
   inic_rast = 2*rand; inic_ugao = 2*pi*rand;
   x = inic_rast*cos(inic_ugao);
   y = inic_rast*sin(inic_ugao);
   ugao_brzine = 2*pi*rand;
   cestice(i,:) = [x, y, ugao_brzine];
end

Ts = 1;
v = 0.5;
std_ro_donja = 0.3; std_ro_gornja = 0.6;
std_teta = pi/36;
najbolje_cestice = [];

%% Česticni filtar
verov_promene_pravca = zeros(N,1);
for t = 1:length(podaci)
    %tezine = 1/N*ones(N,1);
    trenutno_ro = ro_izmereno(t);
    trenutno_teta = teta_izmereno(t);

    for i = 1:N
       p = rand;
       if(p<=verov_promene_pravca(i))
           delta_pravac = -pi/6 + pi/3*rand;
           cestice(i,3) = cestice(i,3) + delta_pravac;  
           verov_promene_pravca(i) = 0.2;
           if cestice(i,3)<0
                  cestice(i,3) = cestice(i,3)+2*pi;
           end
       else
           verov_promene_pravca(i) = verov_promene_pravca(i) + 0.2;
       end
       
       cestice(i,1) = cestice(i,1) + v*cos(cestice(i,3))*Ts + normrnd(0,0.45);
       cestice(i,2) = cestice(i,2) + v*sin(cestice(i,3))*Ts + normrnd(0,0.45);
    end

    for i = 1:N
        teta_estimirano = atan2(cestice(i,2),cestice(i,1)) + laprnd(1,1,0,std_teta);
        std_ro = sum_ro(cestice(i,3), teta_estimirano, std_ro_donja, std_ro_gornja);
        ro_estimirano = sqrt(cestice(i,1)^2+cestice(i,2)^2) + std_ro*randn;
        
        tezina_poteg = exp(-(ro_estimirano - trenutno_ro)^2/(2*std_ro^2))/std_ro/sqrt(pi*2);
        tezina_teta = exp(-abs(teta_estimirano - trenutno_teta)/std_teta)/2/std_teta;
        tezine(i) = tezina_poteg*tezina_teta;
    end
    
    tezine = tezine/sum(tezine);
    %indeksi_reuzorkovani = randsample(N, N, true, tezine);
    [~, indeksi] = sort(tezine,'descend');   
    najbolje_cestice = cestice(indeksi(1:broj_prikazanih_cestica),:);
    najbolje_tezine = tezine(indeksi(1:broj_prikazanih_cestica));
    
    %cestice = cestice(indeksi_reuzorkovani,:);
    
    %verov_promene_pravca = verov_promene_pravca(indeksi_reuzorkovani,:);
    normalizovane_tezine = najbolje_tezine / max(najbolje_tezine);
    
    if t==1 || t==2 
        figure;
        hold on;
        scatter(najbolje_cestice(:,1), najbolje_cestice(:,2), 20 * normalizovane_tezine, 'b', 'filled'); % Najbolje čestice plavim tačkama većim proporcijalno težinama
        plot(trenutno_ro*cos(trenutno_teta), trenutno_ro*sin(trenutno_teta), 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Merenje crvenim krstićem
        hold off;
        title('Prikaz najboljih čestica i merenja za t=2, bez reuzorkovanja');
        xlabel('X koordinata');
        ylabel('Y koordinata');
        legend('Najbolje čestice', 'Merenje');
        if t == 2
        xlim([- 1 1.5])
        ylim([0.5 3])
        else
        end
        grid on
    end

end
sqrt((cestice(1,1))^2+(cestice(1,2))^2);
a = mean(sqrt((cestice(:,1)).^2+(cestice(:,2)).^2));
b = mean(atan(cestice(:,2)./cestice(:,1)));
c = mean(sqrt((najbolje_cestice(:,1)).^2+(najbolje_cestice(:,2)).^2));
d = mean(atan(najbolje_cestice(:,2)./najbolje_cestice(:,1)));


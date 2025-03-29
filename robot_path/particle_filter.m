clear
close all
clc

%% Inicijalizacija

podaci = csvread('observations.csv');
ro_izmereno = podaci(:,1);      % Polarne koordinate - poteg
teta_izmereno = podaci(:,2);    % Polarne koordinate - polarni ugao

broj_prikazanih_cestica = 100;

N = 2000;             % Broj čestica
cestice = zeros(N,3); % Stanje ima 3 komponente, x i y koordinata 
                      % i ugao koji zaklapa vektor brzine čestice sa x-osom 
for i = 1:N
   x = -2 + 4*rand;
   y = -2 + 4*rand;
   ugao_brzine = 2*pi*rand;
   cestice(i,:) = [x, y, ugao_brzine];  % Inicijalizacija
end

Ts = 1;     % Perioda odabiranja
v = 0.5;    % Intenzitet vektora brzine
std_ro_donja = 0.3; std_ro_gornja = 0.6;
std_teta = pi/36;
najbolje_cestice = [];

% Nove promenljive za čuvanje srednjih vrednosti i standardnih devijacija za ro i teta
srednji_ro = zeros(length(podaci), 1);
srednji_teta = zeros(length(podaci), 1);
std_ro_vekt = zeros(length(podaci), 1);
std_teta_vekt = zeros(length(podaci), 1);

% Promenljive za čuvanje putanja
prava_putanja_x = zeros(length(podaci), 1);
prava_putanja_y = zeros(length(podaci), 1);
estimacija_putanja_x = zeros(length(podaci), 1);
estimacija_putanja_y = zeros(length(podaci), 1);

%% Česticni filtar i prikaz prvog i drugog sample-a
verov_promene_pravca = zeros(N,1);
for t = 1:length(podaci)
    tezine = 1/N*ones(N,1);
    trenutno_ro = ro_izmereno(t);
    trenutno_teta = teta_izmereno(t);

    for i = 1:N
       p = rand;
       if(p <= verov_promene_pravca(i))
           delta_pravac = -pi/6 + pi/3*rand;
           cestice(i,3) = cestice(i,3) + delta_pravac;  
           verov_promene_pravca(i) = 0.2;
           if cestice(i,3) < 0
                  cestice(i,3) = cestice(i,3)+2*pi;
           end
       else
           verov_promene_pravca(i) = verov_promene_pravca(i) + 0.2;
       end
       
       cestice(i,1) = cestice(i,1) + v*cos(cestice(i,3))*Ts + normrnd(0,0.3);
       cestice(i,2) = cestice(i,2) + v*sin(cestice(i,3))*Ts + normrnd(0,0.3);
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
    indeksi_reuzorkovani = randsample(N, N, true, tezine);
    [~, indeksi] = sort(tezine,'descend');   
    najbolje_cestice = cestice(indeksi(1:broj_prikazanih_cestica),:);
    najbolje_tezine = tezine(indeksi(1:broj_prikazanih_cestica));
    cestice = cestice(indeksi_reuzorkovani,:);
    
    verov_promene_pravca = verov_promene_pravca(indeksi_reuzorkovani,:);

    % Čuvanje srednjih vrednosti i standardnih devijacija za ro i teta
    ro_vrednosti = sqrt(najbolje_cestice(:,1).^2 + najbolje_cestice(:,2).^2);
    teta_vrednosti = atan2(najbolje_cestice(:,2), najbolje_cestice(:,1));
    srednji_ro(t) = mean(ro_vrednosti);
    srednji_teta(t) = mean(teta_vrednosti);
    std_ro_vekt(t) = std(ro_vrednosti);
    std_teta_vekt(t) = std(teta_vrednosti);

    % Čuvanje pravih putanja
    prava_putanja_x(t) = trenutno_ro * cos(trenutno_teta);
    prava_putanja_y(t) = trenutno_ro * sin(trenutno_teta);
    estimacija_putanja_x(t) = mean(najbolje_cestice(:,1));
    estimacija_putanja_y(t) = mean(najbolje_cestice(:,2));

    if t==1 || t==2     % Prvi i drugi sample (vremenski trenutak)
        figure;
        grid on;
        hold on;
        scatter(najbolje_cestice(:,1), najbolje_cestice(:,2), 200 * najbolje_tezine, 'b', 'filled'); % Najbolje čestice plavim tačkama većim proporcijalno težinama
        plot(trenutno_ro*cos(trenutno_teta), trenutno_ro*sin(trenutno_teta), 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Merenje crvenim krstićem
        hold off;
        title(['Prikaz najboljih čestica i merenja za t=', num2str(t)]);
        xlabel('X koordinata');
        ylabel('Y koordinata');
        legend('Najbolje čestice', 'Merenje');
        if t == 2
        xlim([- 1 1.5])
        ylim([0.5 3])
        else
        end
    end

end

%% Prikazi estimacija za poteg, polarni ugao i putanju robota

% Grafik za ro
figure;
grid on;
hold on;
plot(1:length(podaci), ro_izmereno, 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Merenja crvenim krstićima
plot(1:length(podaci), srednji_ro, 'bo-', 'LineWidth', 2); % Srednje vrednosti plavim kružićima i linijom

% Dodavanje intervala poverenja
for t = 1:length(podaci)
    plot([t t], [srednji_ro(t) - 2*std_ro_vekt(t), srednji_ro(t) + 2*std_ro_vekt(t)], 'g', 'LineWidth', 2);
end

hold off;
grid on;
title('Poređenje merenja i srednjih vrednosti za ro');
xlabel('Vreme (sekunde)');
ylabel('ro');
legend('Merenje', 'Estimacija', 'Interval poverenja (2σ)', 'Location','best');

% Grafik za teta
figure;
hold on;
plot(1:length(podaci), 180*teta_izmereno/pi, 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Merenja crvenim krstićima
plot(1:length(podaci), 180*srednji_teta/pi, 'bo-', 'LineWidth', 2); % Srednje vrednosti plavim kružićima i linijom

% Dodavanje intervala poverenja
for t = 1:length(podaci)
    plot([t t], 180*[srednji_teta(t) - 2*std_teta_vekt(t), srednji_teta(t) + 2*std_teta_vekt(t)]/pi, 'g', 'LineWidth', 2);
end

hold off;
grid on;
title('Poređenje merenja i srednjih vrednosti za teta');
xlabel('Vreme (sekunde)');
ylabel('teta (stepeni)');
legend('Merenje', 'Estimacija', 'Interval poverenja (2σ)', 'Location','best');

% Grafik za putanje
figure;
hold on;
plot(prava_putanja_x, prava_putanja_y, 'r-', 'LineWidth', 2); % Prava putanja crvenom linijom
plot(estimacija_putanja_x, estimacija_putanja_y, 'b-', 'LineWidth', 2); % Usrednjena putanja plavom linijom
hold off;
title('Poređenje prave putanje i procenjene putanje');
xlabel('X koordinata');
ylabel('Y koordinata');
grid on;
legend('Prava putanja', 'Procenjena putanja', 'Location','best');

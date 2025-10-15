clear;clc;close all;

%% Importo datos
data = readmatrix("EDWARD (6).TXT");
ECG = data(:,1);    % ECG [mV]
Resp = data(:,2);   % Respiración 

N = length(data);   % # de datos
fs = 500;           % f d muestreo
f = (0:N-1)*(fs/N); % vector de frecuencias [Hz]
t = (0:N-1)*(N/fs); % vector de tiempo [s]

%% Espectro de f
espec_ECG = fft(ECG);   % Transformada rápida de Fourier

magn_ECG = abs(espec_ECG/N);    % Normalización de magnitud de espectro de ECG, respecto a cantidad de datos y hago valor absoluto
magn_ECG(1) = magn_ECG(1)/2;    % Ajuste para compontenente DC, pa q no salga duplicado

%% Selección de 95% de componentes con > energia

energia = magn_ECG.^2;                  % Energia d cada componente de f d ECG
umbral95 = quantile(energia, 0.95);     % Umbral del 95% de energia
magn_ECG_95 = find(energia >= umbral95);% Magnitud de espectro de ECG con los componentes con 95% de energia

%% Reconstrucción con 95% de componentes con > energia
espec_ECG_reduced = zeros(size(espec_ECG));                  % Inicializo vecto con 0's, del tamaño q quiero
espec_ECG_reduced(magn_ECG_95) = espec_ECG(magn_ECG_95);     % Lleno vector con componentes q suman el 95% de energia, el resto de valores siguen en 0
ECG_recons95 = ifft(espec_ECG_reduced, 'symmetric');         % Reconstrucción

error95 = norm(ECG - ECG_recons95);

%% Graficos
figure;

subplot(3,1,1);         % Grafica de ECG original
plot(t, ECG);
title("ECG Original");
xlabel("Tiempo [s]")
ylabel("[mV]")

subplot(3,1,2);         % Grafica de Espectro de Magnitud de ECG Orignal
plot(f(1:N/2), magn_ECG(1:N/2),'b');
hold on;
plot(f(1:N/2), espec_ECG_reduced(1:N/2),'r');
legend("Original","95 %");
title("Espectro de Magnitud de ECG");
xlabel("Frecuencia [Hz]");
ylabel("Magnitud");

subplot(3,1,3);         % Grafica de ECG reconstruido
plot(t, ECG_recons95);
xlabel("Tiempo [s]")
ylabel("[mV]")

disp(['Error entre la señal original y la reconstruida: ', num2str(error95)]);


% Configuración del test de Montecarlo
numSimulaciones = 10000;
numDatos = length(ECG);
fs = 500; % Frecuencia de muestreo

% Diferencia observada en los datos reales
diffObservada = mean(Cxy) - mean(Cxy_sin10);

% Vector para almacenar las diferencias simuladas
diferenciasSimuladas = zeros(numSimulaciones, 1);

for i = 1:numSimulaciones
    % Generar señales aleatorias (por ejemplo, ruido blanco)
    ECG_sim = randn(numDatos, 1);
    Resp_sim = randn(numDatos, 1);
    
    % Calcular coherencia original
    [Cxy_sim, ~] = mscohere(ECG_sim, Resp_sim, [], [], [], fs);
    
    % Eliminar 10% de datos aleatoriamente
    indices_eliminar = randperm(numDatos, round(0.1*numDatos));
    ECG_sim_sin10 = ECG_sim;
    ECG_sim_sin10(indices_eliminar) = 0;
    Resp_sim_sin10 = Resp_sim;
    Resp_sim_sin10(indices_eliminar) = 0;
    
    % Calcular coherencia después de eliminar datos
    [Cxy_sim_sin10, ~] = mscohere(ECG_sim_sin10, Resp_sim_sin10, [], [], [], fs);
    
    % Calcular diferencia en coherencia
    diferenciasSimuladas(i) = mean(Cxy_sim) - mean(Cxy_sim_sin10);
end

% Calcular el p-valor
p_valor = sum(abs(diferenciasSimuladas) >= abs(diffObservada)) / numSimulaciones;

if p_valor < 0.05
    disp('La diferencia en coherencia es significativa.');
else
    disp('La diferencia en coherencia no es significativa.');
end



%% =========================================================================
% TEST ROBUSTEZZA: RUMORE DI MISURA CON DISTURBO A GRADINO
%% =========================================================================

%% 1. PARAMETRI RUMORE

noise_std = 0.001;      % 1 mm
noise_bandwidth = 100;  % Hz
Ts_simulation = 0.01;   % 10 ms

fprintf('\n========================================\n');
fprintf('TEST ROBUSTEZZA AL RUMORE DI MISURA\n');
fprintf('========================================\n');
fprintf('Rumore std: %.1f mm\n', noise_std*1000);
fprintf('Banda rumore: %d Hz\n', noise_bandwidth);
fprintf('Passo simulazione: %.0f ms\n', Ts_simulation*1000);
fprintf('Disturbo: GRADINO (step) di %.2f m\n', 0.05);
fprintf('========================================\n\n');

%% 2. GENERAZIONE RUMORE

t = 0:Ts_simulation:15;
bump_height = 0.05;  % Altezza gradino [m]

% Genera rumore bianco filtrato (INVARIATO)
rng(42);
Fs = 1/Ts_simulation;
Wn = min(0.9, noise_bandwidth / (Fs/2));

fprintf('Fs = %.1f Hz, Wn = %.4f\n', Fs, Wn);

noise_raw = noise_std * randn(length(t), 1);
[b_filt, a_filt] = butter(2, Wn, 'low');
noise_filtered = filter(b_filt, a_filt, noise_raw);
noise_filtered = noise_filtered * (noise_std / std(noise_filtered));

fprintf('Rumore generato: std = %.3f mm\n\n', std(noise_filtered)*1000);

%% 3. DEFINIZIONE SISTEMI

% Sistema SISO per analisi: Fa → zs-zu (deflessione)
P_def = ss(A, B(:,1), C(2,:), D(2,1));
P_def.InputName = 'Fa';
P_def.OutputName = 'zs-zu';

% Sistema completo per simulazione
sys_full = ss(A, B, C, D);
sys_full.InputName = {'Fa', 'zr'};
sys_full.OutputName = {'zs_ddot', 'zs-zu'};

% Controllore K già sintetizzato
fprintf('Controllore: ordine %d\n', order(K));
fprintf('Loop gain a 1 rad/s: %.2f dB\n', 20*log10(abs(evalfr(P_def*K, 1j))));

%% 4. CALCOLO FUNZIONI DI TRASFERIMENTO

% Loop gain
L = P_def * K;

% Sensibilità
S = feedback(1, L);

% Sensibilità complementare
T = feedback(L, 1);

% Verifica
S_peak = getPeakGain(S);
T_peak = getPeakGain(T);

fprintf('\n--- Proprietà Loop ---\n');
fprintf('Picco di S: %.2f (%.1f dB)\n', S_peak, 20*log10(S_peak));
fprintf('Picco di T: %.2f (%.1f dB)\n', T_peak, 20*log10(T_peak));
fprintf('Norma ∞ di T: %.2f\n', norm(T, inf));
fprintf('→ Amplificazione teorica rumore: %.2fx\n\n', norm(T, inf));

%% 5. SIMULAZIONE SENZA RUMORE - DISTURBO A GRADINO

fprintf('Simulazione senza rumore (disturbo a GRADINO)...\n');

% Chiusura loop: feedback su ingresso 1 (Fa), uscita 2 (zs-zu), negativo
CL_full = feedback(sys_full, K, 1, 2, -1);

% *** MODIFICA CHIAVE: GRADINO invece di IMPULSO ***
% Crea segnale di gradino per ingresso zr (canale 2)
zr_step = bump_height * ones(length(t), 1);  % Gradino permanente da t=0
u_step = [zeros(length(t), 1), zr_step];     % [Fa=0, zr=gradino]

[y_clean, t_out] = lsim(CL_full, u_step, t);

zs_ddot_clean = y_clean(:,1);
zs_zu_clean = y_clean(:,2);

fprintf('✓ Simulazione pulita completata\n');

%% 6. SIMULAZIONE CON RUMORE - METODO SUPERPOSIZIONE

fprintf('Simulazione con rumore (GRADINO + rumore misura)...\n');

% Il rumore sulla misura di zs-zu si propaga attraverso il controllore
% Il rumore genera: Fa_noise = -K * noise
% Questo si propaga attraverso il sistema

% Calcola forza dovuta al rumore
[Fa_from_noise, ~] = lsim(K, -noise_filtered, t);

% Propaga questa forza attraverso il sistema (anello aperto)
% Crea sistema da Fa alle uscite
sys_Fa_to_outputs = ss(A, B(:,1), C, D(:,1));

[y_from_noise, ~] = lsim(sys_Fa_to_outputs, Fa_from_noise, t);

% Risposta totale = risposta gradino + effetto rumore
zs_ddot_noisy = zs_ddot_clean + y_from_noise(:,1);
zs_zu_noisy = zs_zu_clean + y_from_noise(:,2);

% Forza totale
% Calcola Fa del caso pulito
[Fa_clean, ~] = lsim(K, -zs_zu_clean, t);

% Forza totale
Fa_noisy = Fa_clean + Fa_from_noise;

fprintf('✓ Simulazione con rumore completata\n\n');

%% 7. ANALISI COMPARATIVA

fprintf('========================================\n');
fprintf('ANALISI IMPATTO DEL RUMORE\n');
fprintf('========================================\n\n');

% --- Accelerazione ---
fprintf('--- ACCELERAZIONE ---\n');

picco_clean = max(abs(zs_ddot_clean));
rms_clean = rms(zs_ddot_clean);

picco_noisy = max(abs(zs_ddot_noisy));
rms_noisy = rms(zs_ddot_noisy);

noise_in_acc = y_from_noise(:,1);
rms_noise_acc = rms(noise_in_acc);

fprintf('Senza rumore:\n');
fprintf('  Picco: %.4f m/s²\n', picco_clean);
fprintf('  RMS: %.6f m/s²\n', rms_clean);
fprintf('\nCon rumore:\n');
fprintf('  Picco: %.4f m/s²\n', picco_noisy);
fprintf('  RMS: %.6f m/s²\n', rms_noisy);
fprintf('\nEffetto del rumore:\n');
fprintf('  RMS rumore propagato: %.6f m/s²\n', rms_noise_acc);
fprintf('  Amplificazione: %.2fx\n', rms_noise_acc/noise_std);

% --- Deflessione ---
fprintf('\n--- DEFLESSIONE ---\n');

def_max_clean = max(abs(zs_zu_clean));
def_max_noisy = max(abs(zs_zu_noisy));
rms_noise_def = rms(y_from_noise(:,2));

fprintf('Senza rumore: %.2f cm\n', def_max_clean*100);
fprintf('Con rumore: %.2f cm\n', def_max_noisy*100);
fprintf('RMS rumore su deflessione: %.3f mm\n', rms_noise_def*1000);

% --- Forza Attuatore ---
fprintf('\n--- SFORZO DI CONTROLLO ---\n');

rms_Fa_clean = rms(Fa_clean);
rms_Fa_noisy = rms(Fa_noisy);
picco_Fa_clean = max(abs(Fa_clean));
picco_Fa_noisy = max(abs(Fa_noisy));

fprintf('Senza rumore:\n');
fprintf('  RMS: %.2f N\n', rms_Fa_clean);
fprintf('  Picco: %.2f N\n', picco_Fa_clean);
fprintf('\nCon rumore:\n');
fprintf('  RMS: %.2f N\n', rms_Fa_noisy);
fprintf('  Picco: %.2f N\n', picco_Fa_noisy);
fprintf('  Incremento RMS: %+.1f%%\n', (rms_Fa_noisy/rms_Fa_clean-1)*100);
fprintf('  Incremento picco: %+.1f%%\n', (picco_Fa_noisy/picco_Fa_clean-1)*100);

%% 8. VALUTAZIONE ROBUSTEZZA

fprintf('\n========================================\n');
fprintf('VALUTAZIONE ROBUSTEZZA\n');
fprintf('========================================\n\n');

amplification = rms_noise_acc / noise_std;
control_increase = (rms_Fa_noisy/rms_Fa_clean - 1) * 100;
acc_degradation = (picco_noisy - picco_clean) / picco_clean * 100;

fprintf('METRICHE:\n');
fprintf('  Amplificazione rumore: %.2fx', amplification);
if amplification < 10
    fprintf(' ✓✓ Eccellente\n');
elseif amplification < 50
    fprintf(' ✓ Buono\n');
elseif amplification < 100
    fprintf(' → Accettabile\n');
else
    fprintf(' ⚠ Critico\n');
end

fprintf('  Incremento sforzo: %+.1f%%', control_increase);
if control_increase < 20
    fprintf(' ✓ Contenuto\n');
elseif control_increase < 50
    fprintf(' → Moderato\n');
else
    fprintf(' ⚠ Significativo\n');
end

fprintf('  Degrado prestazioni: %+.1f%%\n', acc_degradation);

fprintf('\n--- GIUDIZIO COMPLESSIVO ---\n');

if amplification < 50 && control_increase < 50
    fprintf('✓✓ Sistema ROBUSTO al rumore di misura\n');
    fprintf('   Non necessari filtri aggiuntivi\n');
elseif amplification < 100
    fprintf('✓ Sistema ACCETTABILE\n');
    fprintf('   Considerare filtro solo se prestazioni critiche\n');
else
    fprintf('⚠ Amplificazione ECCESSIVA\n');
    fprintf('   RACCOMANDATO: filtro passa-basso a ~10-20 Hz\n');
end

fprintf('========================================\n\n');

%% 9. ANALISI FREQUENZIALE

fprintf('--- CONFRONTO TEORICO vs SIMULATO ---\n');

% Amplificazione teorica (dalla norma di T)
T_norm_inf = norm(T, inf);

% Amplificazione misurata
amp_measured = amplification;

fprintf('Amplificazione teorica (||T||_∞): %.2f\n', T_norm_inf);
fprintf('Amplificazione misurata: %.2f\n', amp_measured);
fprintf('Rapporto: %.2f%%\n\n', amp_measured/T_norm_inf*100);

fig1 = figure(10);
clf
set(fig1, 'Position', [100, 100, 1200, 600]); % Dimensione figura ottimizzata

% Parametri specifiche
ta1_spec = 3;           % [s]
S_spec = 10;            % [%]
def_max_spec = 0.08;    % [m]
picco_acc_spec = 3;     % [m/s²]

%% CALCOLO METRICHE PER VERIFICA SPECIFICHE

% (1) Errore a regime accelerazione
acc_regime = mean(zs_ddot_noisy(end-100:end));
spec1_ok = abs(acc_regime) < 1e-3;

% (2) Deflessione massima
max_def = max(abs(zs_zu_noisy));
spec2_ok = max_def < def_max_spec;

% (3) Tempo di assestamento
acc_final = mean(zs_ddot_noisy(end-100:end));
toll_1perc = 0.01 * abs(max(zs_ddot_noisy) - acc_final);
if toll_1perc < 1e-6
    toll_1perc = 0.01;
end

upper_band = acc_final + toll_1perc;
lower_band = acc_final - toll_1perc;
out_of_band = (zs_ddot_noisy > upper_band) | (zs_ddot_noisy < lower_band);
idx_last = find(out_of_band, 1, 'last');

if isempty(idx_last)
    ta1_manual = 0;
else
    ta1_manual = t(idx_last);
end
spec3_ok = ta1_manual <= ta1_spec;

% (4) Sovraelongazione
[acc_max, idx_max] = max(zs_ddot_noisy);
[acc_min, idx_min] = min(zs_ddot_noisy);

if abs(acc_max) > abs(acc_min)
    acc_picco = acc_max;
else
    acc_picco = acc_min;
end

if abs(acc_picco) > 1e-6
    S_percent_manual = abs(acc_regime / acc_picco) * 100;
else
    S_percent_manual = 0;
end
spec4_ok = S_percent_manual <= S_spec;

% (5) Picco accelerazione
picco_acc = max(abs(zs_ddot_noisy));
spec5_ok = picco_acc < picco_acc_spec;

% Conteggio specifiche soddisfatte
n_specs_ok = sum([spec1_ok, spec2_ok, spec3_ok, spec4_ok, spec5_ok]);

%% SUBPLOT 1: Accelerazione (più grande, occupa 2/3 della larghezza)
subplot(1,3,[1 2])
plot(t, zs_ddot_noisy, 'r', 'LineWidth', 1.8)
hold on
yline(3, 'k--', 'LineWidth', 1.2)
yline(-3, 'k--', 'LineWidth', 1.2)
yline(acc_final, 'b--', 'LineWidth', 1, 'Color', [0 0.5 0.8])

% Linea verticale del tempo di assestamento con etichetta
if ~isempty(idx_last) && ta1_manual < max(t)
    xline(ta1_manual, 'g--', sprintf('t_a = %.2f s', ta1_manual), ...
        'LineWidth', 2, 'FontSize', 11, 'FontWeight', 'bold', ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom')
end

grid on
xlabel('Tempo [s]', 'FontSize', 12)
ylabel('Accelerazione [m/s²]', 'FontSize', 12)
title('Accelerazione Carrozzeria', 'FontSize', 13, 'FontWeight', 'bold')
xlim([0 15])
ylim([-4 4])
set(gca, 'FontSize', 11)

%% SUBPLOT 2: Deflessione (più grande, occupa 2/3 della larghezza)
subplot(1,3,[1 2])
% Creiamo un secondo asse Y a destra
yyaxis right
plot(t, zs_zu_noisy*100, 'b', 'LineWidth', 1.8)
hold on
yline(8, 'k--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
yline(-8, 'k--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
ylabel('Deflessione [cm]', 'FontSize', 12)
ylim([-10 10])
set(gca, 'YColor', 'b')

% Torniamo all'asse sinistro
yyaxis left
set(gca, 'YColor', 'r')

% In alternativa, usiamo layout separato:
% Rimuoviamo il codice sopra e usiamo questo:

%% LAYOUT CORRETTO: DUE SUBPLOT SEPARATI

% Ridefiniamo la figura
clf
set(fig1, 'Position', [100, 100, 1400, 500]);

% SUBPLOT 1: Accelerazione
subplot(1,3,1)
plot(t, zs_ddot_noisy, 'r', 'LineWidth', 1.8)
hold on
yline(3, 'k--', 'LineWidth', 1.2)
yline(-3, 'k--', 'LineWidth', 1.2)

% Linea verticale del tempo di assestamento con valore
if ~isempty(idx_last) && ta1_manual < max(t)
    xline(ta1_manual, 'g--', sprintf('Ta1=%.2f s', ta1_manual), ...
        'LineWidth', 2.5, 'FontSize', 12, 'FontWeight', 'bold', ...
        'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom')
end

grid on
xlabel('Tempo [s]', 'FontSize', 12)
ylabel('Accelerazione [m/s²]', 'FontSize', 12)
xlim([0 15])
ylim([-4 4])
set(gca, 'FontSize', 11)

% SUBPLOT 2: Deflessione
subplot(1,3,2)
plot(t, zs_zu_noisy*100, 'b', 'LineWidth', 1.8)
hold on
yline(8, 'k--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
yline(-8, 'k--', 'LineWidth', 1.2, 'Color', [0.5 0.5 0.5]);
grid on
xlabel('Tempo [s]', 'FontSize', 12)
ylabel('Deflessione [cm]', 'FontSize', 12)
xlim([0 15])
ylim([-10 10])
set(gca, 'FontSize', 11)

%% SUBPLOT 3: TABELLA RIEPILOGATIVA (PIÙ GRANDE E LEGGIBILE)
subplot(1,3,3)
axis off

% Titolo tabella
text(0.5, 0.98, 'VERIFICA SPECIFICHE H∞', ...
    'FontSize', 14, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top');

% Posizione iniziale
y_pos = 0.88;
line_height = 0.155;

% Array delle specifiche
specs_text = {
    sprintf('(1) Errore regime:\n    %.2e m/s²', acc_regime);
    sprintf('(2) Deflessione max:\n    %.2f cm', max_def*100);
    sprintf('(3) Tempo assest.:\n    %.2f s', ta1_manual);
    sprintf('(4) Sovraelongazione:\n    %.1f%%', S_percent_manual);
    sprintf('(5) Picco accel.:\n    %.2f m/s²', picco_acc);
};

specs_target = {
    '= 0';
    '< 8 cm';
    '≤ 3 s';
    '≤ 10%';
    '< 3 m/s²';
};

spec_ok = [spec1_ok; spec2_ok; spec3_ok; spec4_ok; spec5_ok];

% Stampa le specifiche con formato migliorato
for i = 1:length(specs_text)
    % Rettangolo di sfondo alternato
    if mod(i, 2) == 0
        rectangle('Position', [0.02, y_pos-0.07, 0.96, 0.14], ...
            'FaceColor', [0.95 0.95 0.95], 'EdgeColor', 'none');
    end

    % Testo specifica
    text(0.05, y_pos, specs_text{i}, ...
        'FontSize', 10.5, ...
        'VerticalAlignment', 'top', ...
        'FontName', 'Courier New', ...
        'FontWeight', 'normal');

    % Target
    text(0.65, y_pos, specs_target{i}, ...
        'FontSize', 11, ...
        'VerticalAlignment', 'top', ...
        'FontWeight', 'bold', ...
        'Color', [0.2 0.2 0.2]);

    % Simbolo OK/NON OK (più grande)
    if spec_ok(i)
        text(0.88, y_pos-0.03, '✓', ...
            'FontSize', 18, ...
            'Color', [0 0.7 0], ...
            'FontWeight', 'bold', ...
            'VerticalAlignment', 'top');
    else
        text(0.88, y_pos-0.03, '✗', ...
            'FontSize', 18, ...
            'Color', [0.9 0 0], ...
            'FontWeight', 'bold', ...
            'VerticalAlignment', 'top');
    end

    y_pos = y_pos - line_height;
end

% Linea separatore
y_pos = y_pos + 0.04;
line([0.05 0.95], [y_pos y_pos], 'Color', 'k', 'LineWidth', 2);

% Riepilogo finale
y_pos = y_pos - 0.08;
text(0.5, y_pos, sprintf('Specifiche: %d/5', n_specs_ok), ...
    'FontSize', 13, ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');


%% Salvataggio
nome_file_obs = 'Hinf_control_test_completo.eps';
print(fig1, nome_file_obs, '-depsc', '-painters');
fprintf('✓ Grafico salvato come: %s\n\n', nome_file_obs);

%% STAMPA RIEPILOGO A CONSOLE
fprintf('\n========================================\n');
fprintf('RIEPILOGO VERIFICA SPECIFICHE\n');
fprintf('========================================\n\n');

for i = 1:5
    fprintf('%s\n', strrep(specs_text{i}, '\n    ', ': '));
    fprintf('  Obiettivo: %s\n', specs_target{i});
    if spec_ok(i)
        fprintf('  ✓ SODDISFATTA\n\n');
    else
        fprintf('  ✗ NON SODDISFATTA\n\n');
    end
end

fprintf('========================================\n');
fprintf('Risultato: %d/5 specifiche soddisfatte\n', n_specs_ok);
fprintf('========================================\n\n');
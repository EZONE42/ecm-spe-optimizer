%% run_transference.m — Bruce–Vincent 법으로 Li⁺ 전이수 t⁺ 계산 (M4 농도분극 입력)
% 데이터: Zenodo 7428858, T+_EIS_CA_README.txt 기준
%   대칭셀, 55 °C, EIS 0.1 Hz~1 MHz (10 mV), CA ΔV = 50 mV, README는 4 h
%   (실제 파일: Na 4.0 h, K 6.0 h — 2026-10-29 확인. I_ss는 시간과 무관하게 마지막 10분 평균)
%   전극 면적: Li·Na 0.5026 cm², K 2.0106 cm²
% 식: t+ = I_ss(ΔV − I_0·R_i,0) / ( I_0·(ΔV − I_ss·R_i,ss) )
%     I_0 = ΔV/(R_b,0 + R_i,0) (계산값 권장, 측정 첫 점은 이중층 충전으로 과대),  I_ss = 마지막 10분 평균
% 실행: T+_EIS_*, T+_CA_* 파일과 같은 폴더에서 run_transference
clear; clc;
dV = 0.050;  % V
cases = table(["Li";"Na";"K"], ["T+_EIS_1";"T+_EIS_3";"T+_EIS_5"], ["T+_EIS_2";"T+_EIS_4";"T+_EIS_6"], ...
              ["T+_CA_1";"T+_CA_2";"T+_CA_3"], [0.5026;0.5026;2.0106], ...
              'VariableNames', {'cation','eis_before','eis_after','ca','A_cm2'});
rows = [];
for k = 1:height(cases)
    try
        [Rb0, Ri0] = fitSymmetric(loadEIS(cases.eis_before(k)));
        [Rbs, Ris] = fitSymmetric(loadEIS(cases.eis_after(k)));
        [t, I] = loadCA(cases.ca(k));
        I0  = dV / (Rb0 + Ri0);                       % 계산 초기전류 (A)
        Iss = mean(I(t >= t(end) - 600));             % 마지막 600 s 평균 (A)
        tplus = Iss*(dV - I0*Ri0) / (I0*(dV - Iss*Ris));
        I0pk  = max(I(t <= 5));                        % 측정 초기 최대전류(첫 5 s): 장비 안정화로 0.3 s 무렵 최대
        tplus_pk = Iss*(dV - I0pk*Ri0) / (I0pk*(dV - Iss*Ris));
        rows = [rows; {cases.cation(k), Rb0, Ri0, Rbs, Ris, I0, Iss, tplus, I0pk, tplus_pk, ...
                       Ri0*cases.A_cm2/2, Ris*cases.A_cm2/2}]; %#ok<AGROW>
    catch ME
        warning('%s 처리 실패: %s', cases.cation(k), ME.message);
    end
end
if isempty(rows), error('결과 없음: 파일 이름·열 이름을 확인하세요.'); end
R = cell2table(rows, 'VariableNames', {'cation','Rb0_Ohm','Ri0_Ohm','Rb_ss_Ohm','Ri_ss_Ohm', ...
    'I0_calc_A','Iss_A','t_plus_calcI0','I0_peak_A','t_plus_peakI0','ASR_int0_Ohm_cm2_per_side','ASR_int_ss_Ohm_cm2_per_side'});
R.bulk_fraction = R.Rb0_Ohm ./ (R.Rb0_Ohm + R.Ri0_Ohm);   % ΔV 중 벌크에 걸리는 비율
writetable(R, 'transference_results.csv'); disp(R);
% 판정: bulk_fraction < 0.1 이면 ΔV 대부분이 계면에 걸려 t+가 R 추정 오차에 매우 민감함 → "참고값"으로만 보고
if any(R.bulk_fraction < 0.1)
    warning('벌크 비율 10%% 미만 시료: %s → t+는 신뢰구간이 넓은 참고값', strjoin(R.cation(R.bulk_fraction<0.1), ', '));
end
% 참고: 대칭셀의 계면 R_i는 두 계면의 합 → 한 계면 ASR = R_i·A/2

function E = loadEIS(stem)
% 실제 파일 형식(2026-10-29 확인): 1행 빈 줄, 2행 "Re(Z)/Ohm<TAB>-Im(Z)/Ohm", 주파수 열 없음
% 점 순서 = 고주파 -> 저주파
    c = dir(char(stem + "*")); if isempty(c), error('파일 없음: %s', stem); end
    fn = fullfile(c(1).folder, c(1).name);
    M = readmatrix(fn, 'FileType','text', 'Delimiter','\t', 'NumHeaderLines', 2);
    M = M(all(isfinite(M(:,1:2)),2), :);
    E.re = M(:,1); E.mim = M(:,2);                     % Re(Z), -Im(Z)
end

function [t, I] = loadCA(stem)
% 실제 형식: "time/h<TAB>I/mA"  -> 초, 암페어로 변환
    c = dir(char(stem + "*")); if isempty(c), error('파일 없음: %s', stem); end
    fn = fullfile(c(1).folder, c(1).name);
    hdr = lower(strtrim(fileread(fn))); hdr = hdr(1:min(80,end));
    M = readmatrix(fn, 'FileType','text', 'Delimiter','\t', 'NumHeaderLines', 1);
    M = M(all(isfinite(M(:,1:2)),2), :);
    t = M(:,1); I = M(:,2);
    if contains(hdr, 'time/h'), t = t*3600; elseif contains(hdr, 'time/min'), t = t*60; end
    if contains(hdr, 'i/ma'), I = I*1e-3; elseif contains(hdr, 'i/µa') || contains(hdr, 'i/ua'), I = I*1e-6; end
end

function [Rb, Ri] = fitSymmetric(E)
% 주파수가 없으므로 Nyquist 기하로 구함
%   R_b   = 벌크 반원이 끝나는 점(-Im 첫 극소)의 Re
%   R_tot = 계면 반원(눌린 반원)에 원을 맞춰 실수축과 만나는 오른쪽 점 (Kasa 원 피팅)
    re = E.re; y = E.mim; n = numel(re);
    i1 = find(y(2:end-1) <= y(1:end-2) & y(2:end-1) <= y(3:end), 1) + 1;
    if isempty(i1), [~, i1] = min(y); end
    Rb = re(i1);
    [ymax, k] = max(y(i1:end)); ip = i1 + k - 1;
    iv = n;
    for j = ip+1:n-1
        if y(j) <= y(j-1) && y(j) <= y(j+1) && y(j) < 0.2*ymax, iv = j; break; end
    end
    s = (i1+3):(iv-1); x = re(s); z = y(s);
    c = [x z ones(size(x))] \ (x.^2 + z.^2);
    xc = c(1)/2; yc = c(2)/2; r = sqrt(c(3) + xc^2 + yc^2);
    Rtot = xc + sqrt(max(r^2 - yc^2, 0));
    Ri = Rtot - Rb;
end

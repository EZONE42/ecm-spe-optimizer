%% run_eis_anchor.m — KIT 블록공중합체 전해질(BPE) 온도별 EIS 일괄 피팅 (M3 앵커)
% 데이터: Khudyshkina et al., Zenodo 10.5281/zenodo.7428858 (CC BY 4.0), EIS_README.txt 기준
%   - 9개 파일: BPE-{Li,Na,K}TFSI x EO:A {20,15,10}
%   - 측정: Biologic SP-200, 0.5 Hz ~ 1 MHz, 5 -> 85 -> 5 °C (10 °C 간격), 온도마다 3회
%   - cycle 1-27 = 가열(5..85 °C), cycle 28-54 = 냉각(85..5 °C)
%   - 전극 면적 7.85E-05 m^2 = 0.785 cm^2, 막 두께는 시료별(아래 표)
% 출력: eis_anchor_results.csv (시료·온도·방향·반복별 R_b, σ, 피팅 품질), 그림
% 실행: 이 파일과 EIS_*.txt 파일을 같은 폴더에 두고 MATLAB에서 run_eis_anchor
% 필요: Optimization Toolbox(lsqnonlin). 없으면 USE_FIT = false 로 두면 절편법만 사용
clear; clc;
USE_FIT = true;          % false: 그래프 절편법(R_b = -Z'' 최소점의 Z')만 사용
A_cm2   = 7.85e-5 * 1e4; % m^2 -> cm^2 = 0.785

% ---- 시료 표 (EIS_README.txt 그대로) ----
S = table( ...
  ["EIS_1_BPE-LiTFSI_20-1";"EIS_2_BPE-LiTFSI_15-1";"EIS_3_BPE-LiTFSI_10-1"; ...
   "EIS_4_BPE-NaTFSI_20-1";"EIS_5_BPE-NaTFSI_15-1";"EIS_6_BPE-NaTFSI_10-1"; ...
   "EIS_7_BPE-KTFSI_20-1"; "EIS_8_BPE-KTFSI_15-1"; "EIS_9_BPE-KTFSI_10-1"], ...
  ["Li";"Li";"Li";"Na";"Na";"Na";"K";"K";"K"], ...
  [20;15;10;20;15;10;20;15;10], ...
  [95;140;75;70;145;85;105;140;90], ...
  'VariableNames', {'file','cation','EO_ratio','L_um'});

out = [];
for k = 1:height(S)
    fn = findDataFile(S.file(k));
    if fn == ""
        warning('파일 없음: %s (건너뜀)', S.file(k)); continue;
    end
    D = readBiologic(fn);                     % f, Zre, Zim(부호: 실제 허수부), cyc
    cycles = unique(D.cyc(:))';
    for c = cycles
        idx = D.cyc == c;
        Z = D.Zre(idx) + 1j*D.Zim(idx);
        if ~isempty(D.f)
            f = D.f(idx); [f, o] = sort(f, 'descend'); Z = Z(o);
        end
        [T_C, dirn, rep] = cycleToTemp(c);
        Rb_int = interceptRb(Z);              % 빠른 추정 (스파이크 시작점)
        Rb = Rb_int; n_dl = NaN; chi = NaN; flag = 0;
        if USE_FIT && ~isempty(D.f)
            try
                [Rb, n_dl, chi, flag] = fitBlocking(f, Z, Rb_int);
            catch ME
                warning('%s cycle %d 피팅 실패: %s', S.file(k), c, ME.message);
                Rb = Rb_int; flag = -1;
            end
        end
        L_cm  = S.L_um(k) * 1e-4;
        sigma = L_cm / (Rb * A_cm2);          % S/cm  (σ = L / (R·A))
        out = [out; {S.file(k), S.cation(k), S.EO_ratio(k), S.L_um(k), c, T_C, dirn, rep, ...
                     Rb_int, Rb, sigma, n_dl, chi, flag}]; %#ok<AGROW>
    end
end
if isempty(out), error('처리한 스펙트럼이 없습니다. 파일 이름과 위치를 확인하세요.'); end
R = cell2table(out, 'VariableNames', {'file','cation','EO_ratio','L_um','cycle','T_C', ...
    'direction','rep','Rb_intercept_Ohm','Rb_fit_Ohm','sigma_S_cm','n_dl','chi2_nu','exitflag'});
writetable(R, 'eis_anchor_results.csv');
fprintf('저장: eis_anchor_results.csv (%d개 스펙트럼)\n', height(R));

% ---- 온도별 평균(3회 반복) + 가열/냉각 비교 + Arrhenius 그림 ----
G = groupsummary(R, {'file','cation','EO_ratio','T_C','direction'}, {'mean','std'}, 'sigma_S_cm');
writetable(G, 'eis_anchor_sigma_by_T.csv');
figure('Color','w'); hold on;
mk = containers.Map({'heating','cooling'}, {'o-','s--'});
files = unique(G.file);
for i = 1:numel(files)
    for d = ["heating","cooling"]
        g = G(G.file==files(i) & G.direction==d, :);
        if isempty(g), continue; end
        plot(1000./(g.T_C+273.15), log10(g.mean_sigma_S_cm), mk(char(d)), ...
             'DisplayName', sprintf('%s (%s)', strrep(files(i),'_','\_'), d));
    end
end
xlabel('1000/T (K^{-1})'); ylabel('log_{10} \sigma (S/cm)'); grid on; legend('Location','eastoutside');
title('BPE 전해질 \sigma(T): 가열(실선) vs 냉각(점선)');
saveas(gcf, 'eis_anchor_arrhenius.png');

% ---- Arrhenius / VTF 피팅 (Li 시료, 가열·냉각 각각) ----
fitRows = [];
for i = 1:numel(files)
    for d = ["heating","cooling"]
        g = G(G.file==files(i) & G.direction==d, :);
        if height(g) < 4, continue; end
        T = g.T_C + 273.15; sg = g.mean_sigma_S_cm;
        p = polyfit(1./T, log(sg.*T), 1);               % ln(σT) = lnA − Ea/(kB·T)
        Ea_eV = -p(1) * 8.617e-5;
        vtf = @(b,T) b(1).*T.^-0.5.*exp(-b(2)./(T-b(3)));
        try
            b = lsqcurvefit(vtf, [10 1000 200], T, sg, [1e-6 100 100], [1e8 5000 min(T)-5], ...
                            optimoptions('lsqcurvefit','Display','off'));
        catch
            b = [NaN NaN NaN];
        end
        fitRows = [fitRows; {files(i), d, Ea_eV, b(1), b(2), b(3)}]; %#ok<AGROW>
    end
end
F = cell2table(fitRows, 'VariableNames', {'file','direction','Ea_apparent_eV','VTF_A','VTF_B_K','VTF_T0_K'});
writetable(F, 'eis_anchor_temperature_fits.csv'); disp(F);

%% ===================== 보조 함수 =====================
function fn = findDataFile(stem)
    c = dir(char(stem + "*")); fn = "";
    if ~isempty(c), fn = string(fullfile(c(1).folder, c(1).name)); end
end

function D = readBiologic(fn)
% BioLogic EC-Lab 텍스트 내보내기(탭 구분)를 읽음. 열 이름이 다르면 아래 패턴을 고칠 것
    opts = detectImportOptions(fn, 'FileType','text', 'VariableNamingRule','preserve');
    T = readtable(fn, opts);
    v = string(T.Properties.VariableNames);
    col = @(pat) find(contains(lower(v), pat), 1);
    iF = col('freq'); iRe = col('re(z)'); iIm = col('im(z)'); iC = col('cycle');
    if isempty(iRe) || isempty(iIm)
        error('열을 찾지 못함. 열 이름: %s', strjoin(v, ', '));
    end
    D.Zre = T{:, iRe}; Zim = T{:, iIm};
    if startsWith(strtrim(v(iIm)), '-'), Zim = -Zim; end   % "-Im(Z)" 열이면 부호 되돌림
    D.Zim = Zim; n = numel(D.Zre);
    % 같은 데이터셋의 T+_EIS 파일은 주파수 열 없이 Re/-Im만 있음(2026-10-29 확인) → 대비
    if isempty(iF)
        warning('%s: 주파수 열 없음 → 회로 피팅 대신 절편법(R_b)만 사용', fn);
        D.f = [];
    else
        D.f = T{:, iF};
    end
    if isempty(iC)
        if mod(n, 54) ~= 0, error('cycle 열이 없고 점 수 %d가 54의 배수가 아님', n); end
        D.cyc = repelem((1:54)', n/54);                     % README: 54개 스펙트럼을 같은 길이로 가정
    else
        D.cyc = T{:, iC};
    end
end

function [T_C, dirn, rep] = cycleToTemp(c)
% README: 1-27 가열(5..85 °C, 온도당 3회), 28-54 냉각(85..5 °C)
    if c <= 27
        T_C = 5 + 10*floor((c-1)/3); dirn = "heating"; rep = mod(c-1,3) + 1;
    else
        T_C = 85 - 10*floor((c-28)/3); dirn = "cooling"; rep = mod(c-28,3) + 1;
    end
end

function Rb = interceptRb(Z)
% 차단전극 스펙트럼: 벌크 반원 끝 = 위상각 atan(-Z''/Z')이 최소인 점의 Z' (스파이크 시작점)
% (-Z'' 최소로 고르면 저온에서 거대 벌크 반원의 고주파 시작점을 잘못 고름 — 2026-10-29 실데이터로 확인)
% 실제 파일은 cycle마다 114점 = 57점 왕복 스윕이므로 두 반쪽에서 각각 구해 평균
    n = numel(Z); h = floor(n/2);
    if n >= 20 && mod(n,2) == 0, parts = {Z(1:h), Z(h+1:end)}; else, parts = {Z}; end
    rb = zeros(1, numel(parts));
    for k = 1:numel(parts)
        z = parts{k}; [~, i] = min(atan2(-imag(z), real(z))); rb(k) = real(z(i));
    end
    Rb = mean(rb);
end

function [Rb, n_dl, chi2nu, flag] = fitBlocking(f, Z, Rb0)
% 모델: (R_b || CPE_g) + CPE_dl   (차단전극 SS|SPE|SS 표준 회로)
    w = 2*pi*f(:); Z = Z(:);
    cpe = @(Q,n) 1./(Q.*(1j*w).^n);
    mdl = @(p) (p(1).*cpe(p(2),p(3)))./(p(1)+cpe(p(2),p(3))) + cpe(p(4),p(5));
    islog = logical([1 1 0 1 0]);
    toP = @(x) islog.*10.^x + ~islog.*x;
    wt  = 1./abs(Z);
    res = @(x) [real(mdl(toP(x))-Z).*wt; imag(mdl(toP(x))-Z).*wt];
    p0 = [max(Rb0,1) 1e-10 0.9 1e-6 0.8];
    x0 = p0; x0(islog) = log10(p0(islog));
    lb = x0 - 3; ub = x0 + 3; lb(~islog) = [0.5 0.5]; ub(~islog) = [1 1];
    opt = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',2e4);
    [x, resnorm, r, flag] = lsqnonlin(res, x0, lb, ub, opt);
    p = toP(x); Rb = p(1); n_dl = p(5);
    chi2nu = resnorm / max(numel(r) - numel(x), 1);
end

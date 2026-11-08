import glob
import re

# 모든 EIS txt 파일 검색
txt_files = glob.glob("EIS_*.txt")

for file_path in txt_files:
    output_csv = file_path.replace(".txt", "_formatted.csv")

    formatted_rows = []
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            # 숫자 데이터 3개가 있는 행 추출 (cycle, Re(Z), -Im(Z))
            if len(parts) == 3:
                try:
                    cycle = float(parts[0])
                    re_z = float(parts[1])
                    neg_im_z = float(parts[2])

                    # basic_sample.csv 형식: 3개 열 (Freq/Cycle, Re(Z), Im(Z))
                    # Im(Z)는 -Im(Z)의 음수값 적용
                    formatted_rows.append(f"{cycle},{re_z},{-neg_im_z}\n")
                except ValueError:
                    continue

    with open(output_csv, "w", encoding="utf-8") as out_f:
        out_f.writelines(formatted_rows)

    print(f"변환 완료: {file_path} -> {output_csv}")
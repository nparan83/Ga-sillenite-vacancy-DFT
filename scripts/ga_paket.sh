#!/bin/bash
# ga_paket.sh : SilenitB veri deposu paketi (GitHub + Zenodo). Yalniz kopyalama + gzip; hesap YOK (0 c-saat).
# Kullanim (herhangi bir dizinden):
#   bash /arf/scratch/nparan/SilenitB/99_analiz/ga_paket.sh       -> F5 HARIC tum asamalar
#   bash /arf/scratch/nparan/SilenitB/99_analiz/ga_paket.sh f5    -> yalniz 03_site/f5_* (calc/03_site_f5.tar.gz)
# Cikti: /arf/scratch/nparan/repo_silenitB/{calc,scripts,data} + MANIFEST.md5 ; 9 MB ustu arsivler .part_NN parcalarina bolunur
# Alinmayanlar: out/ ve *.save dizinleri (data-file-schema.xml HARIC), dalga fonksiyonu, yuk yogunlugu, cube, UPF, PAW pdos_atm dosyalari.
set -u
SRC=${GA_SRC:-/arf/scratch/nparan/SilenitB}
OUT=${GA_OUT:-/arf/scratch/nparan/repo_silenitB}
PARCA=9m
MOD=${1:-tum}
cd "$SRC" || { echo "### $SRC yok"; exit 1; }
mkdir -p "$OUT"/calc "$OUT"/scripts/jobs "$OUT"/data

topla() {   # $1 = arsiv adi ; $2 = ek prune deseni (bos olabilir) ; $3.. = kok dizinler
  local ad=$1 pr=$2; shift 2
  local tmp=$OUT/tmp_$ad
  rm -rf "$tmp"; mkdir -p "$tmp"
  local P=( -path '*/out' -o -path '*/out_*' -o -name '*.save' )
  [ -n "$pr" ] && P+=( -o -path "$pr" )
  find "$@" \( "${P[@]}" \) -prune -o -type f \
       \( -name '*.in' -o -name '*.out' -o -name '*.out.*' -o -name '*.dat' -o -name '*.log' \
          -o -name '*.pos' -o -name '*.cif' -o -name '*.txt' \) ! -name '*pdos_atm*' -print0 \
    | xargs -0 -r cp --parents -t "$tmp"
  if [ -n "$pr" ]; then
    find "$@" -path "$pr" -prune -o -path '*.save/data-file-schema.xml' -print0 | xargs -0 -r cp --parents -t "$tmp"
  else
    find "$@" -path '*.save/data-file-schema.xml' -print0 | xargs -0 -r cp --parents -t "$tmp"
  fi
  local n; n=$(find "$tmp" -type f | wc -l)
  if [ "$n" -eq 0 ]; then echo "### $ad: dosya yok"; rm -rf "$tmp"; return; fi
  rm -f "$OUT/calc/$ad".tar.gz "$OUT/calc/$ad".tar.gz.part_*
  tar czf "$OUT/calc/$ad.tar.gz" -C "$tmp" .
  local b; b=$(stat -c %s "$OUT/calc/$ad.tar.gz")
  if [ "$b" -gt 9437184 ]; then
    split -b "$PARCA" -d -a 2 "$OUT/calc/$ad.tar.gz" "$OUT/calc/$ad.tar.gz.part_" && rm -f "$OUT/calc/$ad.tar.gz"
  fi
  printf '%-14s %4d dosya  %6.1f MB (gz)\n' "$ad" "$n" "$(echo "$b/1048576" | bc -l)"
  rm -rf "$tmp"
}

if [ "$MOD" = "f5" ]; then
  ls -d 03_site/f5_* >/dev/null 2>&1 || { echo "### 03_site/f5_* yok"; exit 1; }
  for d in 03_site/f5_*; do grep -q 'JOB DONE' "$d/run.out" 2>/dev/null || echo "### UYARI: $d JOB DONE degil"; done
  topla 03_site_f5 "" 03_site/f5_*
else
  for st in 01_yakinsama 02_O40 03_site 04_eos 05_elek 06_ref 07_hse; do
    [ -d "$st" ] || { echo "### $st yok, atlandi"; continue; }
    topla "$st" '03_site/f5_*' "$st"
  done
  cp 99_analiz/*.py "$OUT/scripts/"
  cp job_*.sh "$OUT/scripts/jobs/" 2>/dev/null
  find 99_analiz -maxdepth 1 -type f \( -name '*.log' -o -name '*.dat' -o -name '*.txt' \) -exec cp -t "$OUT/data/" {} +
  cp 00_pseudo/envanter_B.log "$OUT/data/" 2>/dev/null
  ( cd "$OUT/scripts" && md5sum *.py > scripts_md5.txt )
fi

( cd "$OUT" && find . -type f ! -name MANIFEST.md5 -print0 | sort -z | xargs -0 md5sum > MANIFEST.md5 )
echo "--- $OUT"
du -sh "$OUT"/calc "$OUT"/scripts "$OUT"/data 2>/dev/null
echo "dosya sayisi: $(find "$OUT" -type f | wc -l)  |  9 MB ustu: $(find "$OUT" -type f -size +9M | wc -l) (0 olmali)"

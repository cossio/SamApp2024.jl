#ROOT="/Users/cossio/pCloud Drive/data/2023/SAM-Riboswitch"
ROOT="/Users/jfdcd/projects/2022/SAM_Riboswitch/data/2025-03-03/SAMAP"
rsync -avzP --exclude "*.pdf" --exclude "*.log" --exclude "*.svg" --exclude "*.varna" --delete "$ROOT/" "$ROOT-clean/"
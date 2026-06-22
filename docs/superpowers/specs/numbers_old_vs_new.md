# Audit: benchmark counts old (contaminant/Global.* filtered) vs new (Vadim rule)

_Generated 2026-06-21. Old = report_counts.OLD.tsv (target-only, Global.PG.Q/Global.Q gates).
New = report_counts.tsv (Vadim rule: prot_global=Lib.PG.Q.Value, prec_global=Lib.Q.Value,
prot_perrun_avg/prot_complete=PG.Q.Value, prec_min1/min3=Q.Value; no contaminant/target filter)._

Note: old `proteins_tgt` used Global.PG.Q.Value+target; new `prot_global` uses Lib.PG.Q.Value.
Definitions differ, so deltas mix the filter-removal effect with the q-column switch.

| dataset/version | metric | old | new | Δ | Δ% |
|---|---|--:|--:|--:|--:|
| ProteoBench_Module_7 v1_8_1 | proteins (global total) | 10,957 | 11,296 | +339 | +3.1% |
| ProteoBench_Module_7 v1_8_1 | prot per-run avg | 10,091 | 10,192 | +101 | +1.0% |
| ProteoBench_Module_7 v1_8_1 | prot complete | 8,799 | 8,898 | +99 | +1.1% |
| ProteoBench_Module_7 v1_8_1 | precursors min1 | 116,771 | 126,866 | +10,095 | +8.6% |
| ProteoBench_Module_7 v1_8_1 | precursors min3 | 112,396 | 118,078 | +5,682 | +5.1% |
| ProteoBench_Module_7 v1_8_1 | proteins 2pep | 9,645 | 10,189 | +544 | +5.6% |
| ProteoBench_Module_7 v1_8_1 | peptides | 104,540 | 113,274 | +8,734 | +8.4% |
| ProteoBench_Module_7 v2_5_1 | proteins (global total) | 11,405 | 12,025 | +620 | +5.4% |
| ProteoBench_Module_7 v2_5_1 | prot per-run avg | 11,167 | 11,256 | +89 | +0.8% |
| ProteoBench_Module_7 v2_5_1 | prot complete | 10,270 | 10,348 | +78 | +0.8% |
| ProteoBench_Module_7 v2_5_1 | precursors min1 | 130,246 | 163,598 | +33,352 | +25.6% |
| ProteoBench_Module_7 v2_5_1 | precursors min3 | 129,430 | 148,859 | +19,429 | +15.0% |
| ProteoBench_Module_7 v2_5_1 | proteins 2pep | 10,344 | 10,645 | +301 | +2.9% |
| ProteoBench_Module_7 v2_5_1 | peptides | 118,023 | 121,781 | +3,758 | +3.2% |
| ProteoBench_Module_7 v2_5_1_enterprise | proteins (global total) | 11,587 | 12,229 | +642 | +5.5% |
| ProteoBench_Module_7 v2_5_1_enterprise | prot per-run avg | 11,326 | 11,432 | +106 | +0.9% |
| ProteoBench_Module_7 v2_5_1_enterprise | prot complete | 10,403 | 10,492 | +89 | +0.9% |
| ProteoBench_Module_7 v2_5_1_enterprise | precursors min1 | 133,880 | 168,613 | +34,733 | +25.9% |
| ProteoBench_Module_7 v2_5_1_enterprise | precursors min3 | 133,018 | 151,577 | +18,559 | +14.0% |
| ProteoBench_Module_7 v2_5_1_enterprise | proteins 2pep | 10,463 | 10,809 | +346 | +3.3% |
| ProteoBench_Module_7 v2_5_1_enterprise | peptides | 120,900 | 125,801 | +4,901 | +4.1% |
| PXD049412 v1_8_1 | proteins (global total) | 7,796 | 8,105 | +309 | +4.0% |
| PXD049412 v1_8_1 | prot per-run avg | 7,157 | 7,305 | +148 | +2.1% |
| PXD049412 v1_8_1 | prot complete | 6,124 | 6,246 | +122 | +2.0% |
| PXD049412 v1_8_1 | precursors min1 | 41,756 | 46,205 | +4,449 | +10.7% |
| PXD049412 v1_8_1 | precursors min3 | 40,262 | 43,212 | +2,950 | +7.3% |
| PXD049412 v1_8_1 | proteins 2pep | 6,170 | 6,597 | +427 | +6.9% |
| PXD049412 v1_8_1 | peptides | 39,517 | 43,478 | +3,961 | +10.0% |
| PXD049412 v2_5_1 | proteins (global total) | 8,012 | 8,629 | +617 | +7.7% |
| PXD049412 v2_5_1 | prot per-run avg | 7,895 | 8,044 | +149 | +1.9% |
| PXD049412 v2_5_1 | prot complete | 7,089 | 7,193 | +104 | +1.5% |
| PXD049412 v2_5_1 | precursors min1 | 46,243 | 58,208 | +11,965 | +25.9% |
| PXD049412 v2_5_1 | precursors min3 | 45,880 | 52,292 | +6,412 | +14.0% |
| PXD049412 v2_5_1 | proteins 2pep | 6,554 | 6,760 | +206 | +3.1% |
| PXD049412 v2_5_1 | peptides | 43,766 | 45,601 | +1,835 | +4.2% |
| PXD049412 v2_5_1_enterprise | proteins (global total) | 8,252 | 8,750 | +498 | +6.0% |
| PXD049412 v2_5_1_enterprise | prot per-run avg | 8,130 | 8,190 | +60 | +0.7% |
| PXD049412 v2_5_1_enterprise | prot complete | 7,349 | 7,299 | -50 | -0.7% |
| PXD049412 v2_5_1_enterprise | precursors min1 | 47,807 | 59,694 | +11,887 | +24.9% |
| PXD049412 v2_5_1_enterprise | precursors min3 | 47,422 | 53,369 | +5,947 | +12.5% |
| PXD049412 v2_5_1_enterprise | proteins 2pep | 6,656 | 6,908 | +252 | +3.8% |
| PXD049412 v2_5_1_enterprise | peptides | 45,087 | 47,468 | +2,381 | +5.3% |
| PXD062685 v1_8_1 | proteins (global total) | 10,187 | 10,682 | +495 | +4.9% |
| PXD062685 v1_8_1 | prot per-run avg | 10,199 | 10,302 | +103 | +1.0% |
| PXD062685 v1_8_1 | prot complete | 9,601 | 9,700 | +99 | +1.0% |
| PXD062685 v1_8_1 | precursors min1 | 97,149 | 102,528 | +5,379 | +5.5% |
| PXD062685 v1_8_1 | precursors min3 | 96,617 | 101,430 | +4,813 | +5.0% |
| PXD062685 v1_8_1 | proteins 2pep | 8,997 | 9,325 | +328 | +3.6% |
| PXD062685 v1_8_1 | peptides | 87,605 | 92,328 | +4,723 | +5.4% |
| PXD062685 v2_5_1 | proteins (global total) | 11,245 | 11,805 | +560 | +5.0% |
| PXD062685 v2_5_1 | prot per-run avg | 11,231 | 11,344 | +113 | +1.0% |
| PXD062685 v2_5_1 | prot complete | 10,629 | 10,744 | +115 | +1.1% |
| PXD062685 v2_5_1 | precursors min1 | 117,396 | 137,446 | +20,050 | +17.1% |
| PXD062685 v2_5_1 | precursors min3 | 116,394 | 130,035 | +13,641 | +11.7% |
| PXD062685 v2_5_1 | proteins 2pep | 10,081 | 10,348 | +267 | +2.6% |
| PXD062685 v2_5_1 | peptides | 106,698 | 110,280 | +3,582 | +3.4% |
| PXD062685 v2_5_1_enterprise | proteins (global total) | 11,405 | 11,918 | +513 | +4.5% |
| PXD062685 v2_5_1_enterprise | prot per-run avg | 11,328 | 11,400 | +72 | +0.6% |
| PXD062685 v2_5_1_enterprise | prot complete | 10,729 | 10,784 | +55 | +0.5% |
| PXD062685 v2_5_1_enterprise | precursors min1 | 119,726 | 140,620 | +20,894 | +17.5% |
| PXD062685 v2_5_1_enterprise | precursors min3 | 118,761 | 132,892 | +14,131 | +11.9% |
| PXD062685 v2_5_1_enterprise | proteins 2pep | 10,131 | 10,458 | +327 | +3.2% |
| PXD062685 v2_5_1_enterprise | peptides | 108,670 | 113,236 | +4,566 | +4.2% |
| PXD070049 v1_8_1 | proteins (global total) | 9,303 | 9,527 | +224 | +2.4% |
| PXD070049 v1_8_1 | prot per-run avg | 8,788 | 8,877 | +89 | +1.0% |
| PXD070049 v1_8_1 | prot complete | 7,828 | 7,933 | +105 | +1.3% |
| PXD070049 v1_8_1 | precursors min1 | 88,179 | 96,513 | +8,334 | +9.5% |
| PXD070049 v1_8_1 | precursors min3 | 86,662 | 93,469 | +6,807 | +7.9% |
| PXD070049 v1_8_1 | proteins 2pep | 7,770 | 8,259 | +489 | +6.3% |
| PXD070049 v1_8_1 | peptides | 73,267 | 79,985 | +6,718 | +9.2% |
| PXD070049 v2_5_1 | proteins (global total) | 9,592 | 10,234 | +642 | +6.7% |
| PXD070049 v2_5_1 | prot per-run avg | 9,470 | 9,586 | +116 | +1.2% |
| PXD070049 v2_5_1 | prot complete | 8,806 | 8,926 | +120 | +1.4% |
| PXD070049 v2_5_1 | precursors min1 | 102,457 | 128,099 | +25,642 | +25.0% |
| PXD070049 v2_5_1 | precursors min3 | 101,579 | 118,149 | +16,570 | +16.3% |
| PXD070049 v2_5_1 | proteins 2pep | 8,557 | 8,906 | +349 | +4.1% |
| PXD070049 v2_5_1 | peptides | 86,787 | 91,344 | +4,557 | +5.3% |
| PXD070049 v2_5_1_enterprise | proteins (global total) | 9,599 | 10,230 | +631 | +6.6% |
| PXD070049 v2_5_1_enterprise | prot per-run avg | 9,634 | 9,735 | +101 | +1.0% |
| PXD070049 v2_5_1_enterprise | prot complete | 8,938 | 9,011 | +73 | +0.8% |
| PXD070049 v2_5_1_enterprise | precursors min1 | 105,639 | 132,254 | +26,615 | +25.2% |
| PXD070049 v2_5_1_enterprise | precursors min3 | 104,709 | 120,929 | +16,220 | +15.5% |
| PXD070049 v2_5_1_enterprise | proteins 2pep | 8,607 | 8,950 | +343 | +4.0% |
| PXD070049 v2_5_1_enterprise | peptides | 88,987 | 93,777 | +4,790 | +5.4% |

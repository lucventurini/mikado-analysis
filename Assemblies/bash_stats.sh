#bsub -J lessthan1 -q Prod128  -o out_grep -R "span[ptile=4]rusage[mem=2000]"  -n 4 " awk  '{if ($13<1) print $0}' coverage_intron.gff > coverage_less_than_1_intron.txt " ; 
#system ("bsub -w 'ended(\"$jobname\_coveragebed\")' -J $jobname\_Transcript_ids -q $queue  -o out_grep -R \"span[ptile=2]rusage[mem=5000] \" -n 2 \" awk '{if (\$3==\"gene\")print \$9}' $infile1|sed 's/;/\t/g;s/ID=//g'|sort -u > List_transcript_ids.txt\"");

bsub -q Test128 -o out_awk_intron1 -J coverage_less_intron "awk  '{if (\$13<1) print \$0}' coverage_intron.gff > coverage_less_than_1_intron.txt";
bsub -q Test128 -o out_for_mikado -J formikado "awk '{if (\$3==\"mRNA\") print \$9}' Reference.added_introns.gff > Ninth_column_ref.txt";

bsub -q Test128 -o out_awk_exon1 -J coverage_less_exon "awk  '{if (\$13<1) print \$0}' coverage_exon.gff > coverage_less_than_1_exon.txt";

bsub -w 'ended("coverage_less_exon")' -J Lessthan1ids_exons -q Test128 -o out_lessthan1ids "awk '{print \$9}' coverage_less_than_1_exon.txt |sed 's/;/\t/g;s/Parent=//g' |awk '{print \$2}'|sort -u > List_of_transcripts_less_than_1_coverage_exon.txt ";
bsub -w 'ended("Lessthan1ids_exons")' -J comm_compare_exon -q Test128 -o out_comm13 "comm -13 List_of_transcripts_less_than_1_coverage_exon.txt List_transcript_ids.txt > List_of_transcript_ids_100_percent_covered_exon.txt";



bsub -q Test128 -o out_awk_list_transcripts "awk '{if (\$3==\"mRNA\")print \$9}' Reference.added_introns.gff |sed 's/;/\t/g;s/ID=//g'|awk '{print \$1}'|sort -u > List_transcript_ids.txt";

bsub -w 'ended("coverage_less_intron")' -J Lessthan1ids_introns -q Test128 -o out_lessthan1ids "awk '{print \$9}' coverage_less_than_1_intron.txt |sed 's/Parent=//g' |sort -u > List_of_transcripts_less_than_1_coverage_intron.txt ";

bsub -w 'ended("Lessthan1ids_introns")' -J comm_compare_intron -q Test128 -o out_comm13 "comm -13 List_of_transcripts_less_than_1_coverage_intron.txt List_transcript_ids.txt > List_of_transcript_ids_100_percent_covered_intron.txt";

bsub -w 'ended("comm_compare_exon")' -J comm_compare_intron_exon_cov -q Test128 -o out_comm12 " comm -12 List_of_transcript_ids_100_percent_covered_exon.txt List_of_transcript_ids_100_percent_covered_intron.txt > Transcripts_covered_byintron_exon_100percent.txt";

bsub -w 'ended("comm_compare_intron_exon_cov")' -J grepformikado -q Test128 -o out_list_mikado_grep "grep -f Transcripts_covered_byintron_exon_100percent.txt  Ninth_column_ref.txt > formikado_grep.txt ";
bsub -w 'ended("grepformikado")' -J convertformikado -q Test128 -o out_mikadogrepconvert "sed 's/;/\t/g;s/ID=//g;s/Parent=//g' formikado_grep.txt |awk '{print \$1\"\t\"\$2}' > Mikado_reffilter.txt ";
bsub -w 'ended("convertformikado")' -J Ref_filter -q Test128 -o Mikadoreffilter  "source mikado-0.9.3 ;mikado.py util grep Mikado_reffilter.txt TAIR10_GFF3_genes.gff  FilteredREF_AT.gff"
#ID=FBtr0299529;geneID=FBgn0040037
#ID=AT1G01040.2.exon7;Parent=AT1G01040.2

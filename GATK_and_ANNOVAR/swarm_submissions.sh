swarm -g 32 -f scripts/MarkDuplicates.sh --time 48:00:00 --module picard
# 14214483

swarm -g 32 --gres=lscratch:10 --time 48:00:00 --module picard,GATK -f scripts/Index_Realigner.sh
# 14240888

swarm -g 32 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/Realign_Caller.sh
# 14250124 ## gatk had trouble with pseudochromosomal sequences so I will use a targets file

swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/Caller.sh
# 14293079 - cancelled since it looks like chrEBV messed up the indel realignment step before.

swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/Caller.sh
# so actually indel realignment is unnecessary for the HaplotypeCaller. going directly from deduped reads to vcf.
# 14294915 -- cancelled- need to do BSQR.

swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/BSQR.sh
# 14313944

swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/SplitVcf.sh #14395150

swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/FilterVcf.sh #14687477

swarm -g 2 --gres=lscratch:2 --time 48:00:00 --module annovar -f scripts/annovar_run.sh  # ran twice, once for each sample individually and once for a merged vcf. # 14734810


swarm -g 16 --gres=lscratch:10 --time 48:00:00 --module GATK -f scripts/FilterVcf2.sh #


swarm -g 2 --gres=lscratch:2 --time 48:00:00 --module annovar -f scripts/annovar_run2.sh

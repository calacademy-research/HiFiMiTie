#!/bin/bash

# get the CR lengths from the CR excerpt file using the CR length excluding any flanking region length
# get a sorted list of the lengths with counts by uniq and place comments between lengths
# where the jumo is greater the 10 bases.

# from this we will define blocks of lengths that should define different CRs if there are any

# we need at least 20 (for now) to define one of these and if we have a set thst has more than
# 50 recs then will remove the singeltons from both ends

: '  example header
>m64049_220630_141611/179045291/ccs_RC_Val_CR_Met 4706..5988 1283nt Val 66 Met 68 CR 1149
'

BLOCK_SEP=12  # separate when lengths jump this much or more
BLOCK_MIN=3  # need at least this many reads

function make_CR_lengths_file {
   local dir=$1
   [ -z $dir ] && dir="."

   CR_fasta_file=CR_recs_w_flanks.fasta

   # make sure we are getting numbers from header lines stating clearly that that number is for the CR
   awk '/^>.*CR [0-9]+$/{print $NF}' ${dir}/CR_recs_w_flanks.fasta | sort -k2,2r | uniq -c |
   awk -v BLOCK_SEP=$BLOCK_SEP 'NR==1{lst=$2;print; next}$2-lst > BLOCK_SEP{print "# " $2-lst}{print; lst = $2}' > ${dir}/CR_lengths.srt # put comments between close blocks
}

function CR_blocks_from_lengths {
   lengths_file=$1
   awk -v BLOCK_MIN=$BLOCK_MIN '
      BEGIN{ blk = 1; n = 0  }
      /^#/{ blk++; n = 0; next }

      { len[blk][++n] = $2; len_count[blk][n] = $1; blkrecs[blk] += $1 }

      END {
         for (b = 1; b <= blk; b++) {
            nlen = length(len[b])
            if (blkrecs[b] < BLOCK_MIN)
               ignore_block(b)
            else
               cull_and_print_block(b)
         }
      }
      function ignore_block(blkno,     i) {
         printf ("# %d read(s), block, %d too few reads:", blkrecs[blkno], blkno)
         for (i=1; i <= length(len[blkno]); i++)
            printf(" %d of %d", len_count[blkno][i], len[blkno][i])
         printf("\n")
      }
      function cull_and_print_block(blkno,    i, nlen, c_short, c_long, c_short_reads, c_long_reads) {
          nlen = length(len[blkno]); nreads = blkrecs[blkno]
          printf("# %d reads of lengths %d to %d Block %d", nreads, len[blkno][1], len[blkno][nlen], blkno)

          # culling at start at ends if we have a lot of them
          cull_min = nreads; c_short_reads = 0; c_long_reads = 0; c_long = nlen
          if (nreads > 100) cull_min = 1 # cull at least singletons
          if (nreads > 1000) cull_min = 2 # cull singles or doubles

          for (c_short=1; c_short <= nlen && len_count[blkno][c_short] <= cull_min; c_short++) {  # cull till c_short-1 at beginnning
             c_short_reads += len_count[blkno][c_short]
          }
          for (c_long=nlen; c_long >= 1 && len_count[blkno][c_long] <= cull_min; c_long--) {
             c_long_reads += len_count[blkno][c_long]
          }
          if (c_short_reads || c_long_reads) {
             for (i = c_short; i <= c_long; i++) {
                new_read_count += len_count[blkno][i]
             }
             if (new_read_count < 100) { # do not do the culling
                c_short = 1; c_long = nlen; new_read_count = nreads
             }
             else printf(". culling %d reads at start %d at end for block of %d reads lengths %d to %d",
                         c_short_reads, c_long_reads, new_read_count, len[blkno][c_short], len[blkno][c_long])
          }
          dist_in_block = len[blkno][c_long] - len[blkno][c_short]
          printf("\n")
          print len[blkno][c_short], len[blkno][c_long], new_read_count, "reads distance", dist_in_block
      }
   ' $lengths_file
}

function msg {
   echo -e "$@" >&2
}
function set_args {
   dir=$1
   [ -z $dir ] && dir="."

   [ ! -z "$2" ] && BLOCK_SEP=$2
   [ ! -z "$3" ] && BLOCK_MIN=$3

   msg "reading $dir/CR_recs_w_flanks.fasta for CR read blocks with min length difference of $BLOCK_SEP and min number of reads $BLOCK_MIN"

   output_file=$dir/CR_block_read_lengths_max_sep_${BLOCK_SEP}_min_recs_${BLOCK_MIN}.info
   msg "usage CR_length_blocks.sh <dir> <min separator> <min reads>"
   msg "\noutput in $output_file"
}

set_args $@
make_CR_lengths_file $dir
CR_blocks_from_lengths CR_lengths.srt > $output_file

#!/bin/bash

# using the mifi file, show the counts of the trnas to the right of a trna in each of the mito recs
# change to have a minimum for most_freq to use it and also ignore X

# 07Sep2022  change way gh is handled since there can be more than one we need to differntiate or it causes problems

mitfi=$1
[ ! -s "$mitfi" ] && mitfi=mito_hifi_recs.mitfi

show_first=$2
[ -z "$show_first" ] && show_first="F"

function right_neighbor_matrix {
  awk -v show_first=$show_first '
    BEGIN{OFS="\t"; show_first = toupper(show_first); sp = " "}

    function most_freq_neighbor(t) {
      biggest_freq = 0; most_freq = ""
      if (t in right_neighbor) {
         for (n in right_neighbor[t]) {
            if (right_neighbor[t][n] > biggest_freq) {
               biggest_freq = right_neighbor[t][n]
               most_freq = n
            }
         }
      }
      return most_freq
    }

    function prt_order() {
      for(i=1; i<=o; i++)
         printf("%s ", order[i])
      printf("\n")
    }

    function insert_val_before_this_in_order(val, rt_neighbor) {
      num_items=length(order)
      n = 1
      for(i=1; i<=num_items; i++) { # copy items before we see rt_neighbor into new array
         if(order[i] != rt_neighbor)
           nord[n++] = order[i]
         else
            break
      }
      if(order[i] == rt_neighbor) # insert val here
         nord[n++] = val
      for(; i<=num_items; i++) {
         nord[n++] = order[i]
      }
      delete order

      nord_len = length(nord)
      for(n=1; n <= nord_len; n++)
         order[n] = nord[n]
    }

    # skip comment and low quality hits
    { skip = /^#/ || ! ($5 ~ "E") || $7 == "X" }

    # handle goosehairpin and 12S, 16S rrnas
    $8 == "goose_hairpin" {
      skip = 0
      $7 = "gh"
    }
    $8 ~ "^12S" {
      skip = 0
      $7 = "12S"
    }
    $8 ~ "^16S" {
      skip = 0
      $7 = "16S"
    }
    $8 ~ "^OL" {
      skip = 0
      $7 = "OL"
    }
    $7 == "OH" {
      skip = 0
    }

    # skip comment and low quality hits
    skip {next}

    # new record, set a couple vars and get next line
    $1 != lst_rec { cur=$7; lst=""; lst_rec=$1; next }

    {
      right_neighbor[cur][$7]++

      lst = cur
      cur = $7
      trna[cur]++
    }

    END {
      # fill array in the order seen most often
      order[++o] = show_first; delete trna[show_first]
      for(cur=show_first; cur in right_neighbor; cur = mf) {
         mf = most_freq_neighbor(cur)
         if (mf == show_first || mf in seen) {
            break
         }
         order[++o] = mf
         seen[mf]++
         delete trna[mf]
      }
      for (t in trna) { # any we did not hit naturally, see what their neighbor would be
         # order[++o] = t
         mf = most_freq_neighbor(t)
         if(mf!="" && biggest_freq >= 8) { # insert it into the order before the mf value
            if(mf == show_first) # append to end
               order[++o] = t
            else # insert t before mf
               insert_val_before_this_in_order(t, mf)
            o = length(order)
         }
         delete trna[t]
      }

      # print header line in the order seen most often
      for (t1=1; t1 <= o; t1++) {
         printf("\t%s", order[t1])
      }
      printf("\n")

      # now for each one show the counts
      for (t1=1; t1 <= o; t1++) {
         cur = order[t1]
         printf ("%s", cur)
         for (t2=1; t2 <= o; t2++) {
            neigh = order[t2]
            val = int(right_neighbor[cur][neigh])
            val = val==0 ? "." : val
            printf("\t%s", val)
         }
         printf("\n")
      }
    }
  ' <( grep -v -e OH -e goose_hairpin $mitfi)
}

### once we have matrix with the basics and maybe bit more, add the items that can occur more than once  ###
### starting with gh -- we do the right neighbor counts and add the row/cols if they are big enough to   ###
### validate, if no gh in the input mitfi the this will just output as is                                ###

# special item could be gh or OH and perhaps even OL later
function add_special_item {  # reads matrix from stdin
   special_item=$1; [ -z $special_item ] && special_item=gh

   awk -v special_item=$special_item '
      BEGIN { delete valid; delete valid_order }  # delete to make sure valid and valid_order these are treated as arrays

      ###### bevy of functions to manipulate arrays since we need to read all lines before doing the work in the END section ######

      function validate_and_order(mean,   r,v,count) {
         threshold = (mean + 0.5) / 2
         PROCINFO["sorted_in"] = "@val_num_desc"
         for (r in right_neighbor) {  # consider valid if count is 50% or more of mean
            count = right_neighbor[r]
            if (count >= threshold) {
               right_neighbors_ordered_by_count[++v] = r
               valid[r] = count
               right_neighbor_at_pos[itm_pos[r]] = r
            }
         }
         if (length(valid) > 0) {  # make array with valid entries ordered by their predecessor position in descending position order
            for (v in valid) {  # v is the name like "F", vpredpos is int, for first itm it is number of itms, others it is itms pos - 1
               vpredpos = itm_pred_pos[v]
               valid_pred_pos[v] = vpredpos
            }
            PROCINFO["sorted_in"] = "@val_num_desc"
            for (v in valid_pred_pos)
               valid_order[++v_order] = v
         }
         return length(valid)
      }

      function is_valid_pred_pos(m_pos,  v) { # see if m_pos is the value in the the valid_pred_pos array
         for (v in valid_pred_pos) {
            if (valid_pred_pos[v] == m_pos)
               return 1
         }
         return 0
      }


      function arycat(ary, fld1, fld2,       rslt,f) {  # concatenate fields fld1 thru fld2 with tab separator
         rslt = ary[fld1]
         for (f = 2; f <= fld2; f++)
            rslt = rslt "\t" ary[f]
         return rslt
      }
      function join_all(ary) {
         return arycat(ary, 1, length(ary))
      }

      function insert_at(ary, insert_pos, val,        p,new_end,num_to_move) { # ary keys are number from 1 to length(ary)
         if (insert_pos < 0) insert_pos = 1
         if (insert_pos > length(ary))  # this is an append
            append(ary, val)
         else { # move the items from last entry to insert_pos to make room, then can change val of insert_pos
            new_end = length(ary)+1; num_to_move = length(ary) - insert_pos + 1
            for (p = new_end; num_to_move--; p--)
               ary[p] = ary[p-1]
            ary[insert_pos] = val
         }
      }
      function append(ary,     val) {
         ary[ length(ary)+1 ] = val
      }
      function insert_after(ary, insert_pos, val) { insert_at(ary, insert_pos + 1, val)  }

      function insert_into_line_at(line, insert_pos, val,       arry) {
         split(line, arry, "\t")
         insert_at(arry, insert_pos, val)
         return join_all(arry)
      }
      function insert_into_line_after(line, insert_pos, val) { return insert_into_line_at(line, insert_pos + 1, val) }

      function dots(num,    i) {
         dot_str = ""
         for (i=1; i <= num; i++)
            dot_str = dot_str "\t."
         return dot_str
      }

      ##### Start processing lines here ######

      /^--/ || /^#/ || /^$/ {next}

      ## File 1: element (gh, etc) with right neighbor lines grep output handled here ##

      FNR==NR && $7 == special_item { lst_gh = 1; gh_rec = $1; next }
      FNR==NR && lst_gh && $1==gh_rec {  lst_gh = 0; right_neighbor[$7]++ }
      FNR==NR { next }

      ## File 2: matrix file lines handled here ##

      # header line, use it to see order of things
      FNR==1 { # remember header and create 4 arrays for itm order and to get the pos for a given itm, and also succ and pred arrays
         header = $0
         tot_itms = split(substr(header,2), itm_order, "\t")
         for (itm in itm_order)
            itm_pos[ itm_order[itm]  ] = itm
         for (i=1; i <= tot_itms; i++) {
            itm = itm_order[i]
            itm_succ_pos[itm] = (i < tot_itms) ? i+1 : 1; itm_succ_col[itm] = itm_succ_pos[itm] + 1  # itm is in first col so col loc is one greater than pos
            itm_pred_pos[itm] = (i > 1) ? i-1 : tot_itms; itm_pred_col[itm] = itm_pred_pos[itm] + 1  #                ""              ""
         }
      }

      FNR > 1 { # matrix row lines
         matrix[++mlines] = $0
         matrix_itm[mlines] = $1
         matrix_itm_pos[$1] = mlines

         count = (int($2)) ? int($2) : $(FNR+1)
         sum += count
      }

      END {
         mean = sum/tot_itms
         num_valid = validate_and_order(mean) # put those that pass the means test in an array sorted by largest counts. arrays named valid and valid_order created

         for (i = 1; i <= num_valid; i++) {  # insert v in the proper place in the header
            v = valid_order[i]  # valid_order goes from those closest to end to those closest to beginning so we can rely on original col indexes for inserts
            insert_after_pos = itm_pred_col[v]
            header = insert_into_line_after(header, insert_after_pos, special_item)
         }

         print header

         for (m = 1; m <= mlines; m++) {
            if (num_valid == 0) { # no work to do, just print the matrix input
               print matrix[m]
               continue
            }

            newline = ""  # created if we need to add a special item line
            n = split(matrix[m], m_arry, "\t")  # manipulate array, rebuild into line to print at bottom of loop
            m_itm = m_arry[1]

            for (i = 1; i <= num_valid; i++) {  # insert v in the proper place
               v = valid_order[i]
               this_valids_pred = (is_valid_pred_pos(m) && m == itm_pred_pos[v])

               to_insert = "."
               if (this_valids_pred) {  # move the value in the current successor to the inserted successor and replace it with a dot
                  val_to_move = itm_succ_col[m_itm]
                  to_insert = m_arry[val_to_move]
                  m_arry[val_to_move] = "."
               }

               insert_after_pos = itm_pred_col[v]
               insert_after(m_arry, insert_after_pos, to_insert)

               if (this_valids_pred) { # create the line to be inserted for the special item
                  succ_pos = itm_succ_pos[m_itm]; succ_itm = itm_order[succ_pos]
                  dots_before_rn = succ_pos > 1 ? succ_pos : 0
                  newline = special_item dots(dots_before_rn) "\t" valid[succ_itm] dots(tot_itms - dots_before_rn + length(valid) - 1)
               }
            }

            print join_all(m_arry)
            if (newline != "")
               print newline
         }  # end of loop thru matrix lines
      }
   ' <(echo "#so file is not empty" && grep -A1 $special_item $mitfi) -
}

function add_gh {
   add_special_item gh
}

right_neighbor_matrix | add_gh

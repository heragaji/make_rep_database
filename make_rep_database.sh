#! /bin/bash
# ./make_rep_database.sh -i FASTA -d output_dir 
cd `dirname $0`
src_dir=`pwd`/src
PROGNAME=$(basename $0)
VERSION="0.1"
HELP_MSG="'$PROGNAME -h'と指定することでヘルプを見ることができます"

# ヘルプメッセージ
usage() {
  echo "Usage: $PROGNAME -a arg [-b arg] param"
  echo 
  echo "オプション:"
  echo "  -h, --help"
  echo "      --version"
  echo "  -d, --dir <ARG>     output directory"
  echo "  -in, --input <ARG>     input file"
  echo "  -re, --realigner <ARG>     realigner directory"
  echo "  -du, --dump-consensus <ARG>     dump-consensus directory"
  echo "  -top, --top <ARG>     the number of top"
  echo "  -it, --iteration <ARG>     the maximum number of iteration of realigner"
  echo "  -cov, --coverage <ARG>     threshold of depth"
  echo "  -int, --interval <ARG>     threshold of peak intarval"
  echo "  -pe, --peak <ARG>     threshold of peak coverage"
  echo "  -cut, --cut <ARG>     threshold of cut interval"
  echo
  exit 1
}

# オプション解析
for OPT in "$@"
do
  case "$OPT" in
    # ヘルプメッセージ
    '-h'|'--help' )
      usage
      exit 1
      ;;
    # バージョンメッセージ
    '--version' )
      echo $VERSION
      exit 1
      ;;
    # オプション-d、--dir
    '-d'|'--dir' )
      output_dir=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_D="$2"
      shift 2
      ;;
    # オプション-in、--input
    '-in'|'--input' )
      reads_data=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_IN="$2"
      shift 2
      ;;
    # オプション-re、--realigner
    '-re'|'--realigner' )
      realigner_dir=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_RE="$2"
      shift 2
      ;;
    # オプション-du、--dump-consensus
    '-du'|'--dump-consensus' )
      dump_dir=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_DU="$2"
      shift 2
      ;;
    # オプション-top、--top
    '-top'|'--top' )
      top=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_TOP="$2"
      shift 2
      ;;
    # オプション-it、--itration
    '-it'|'--itration' )
      it=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_IT="$2"
      shift 2
      ;;
    # オプション-cov、--coverage
    '-cov'|'--coverage' )
      cov=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_COV="$2"
      shift 2
      ;;
    # オプション-int、--interval
    '-int'|'--interval' )
      interval=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_INT="$2"
      shift 2
      ;;
    # オプション-pe、--peak
    '-pe'|'--peak' )
      pe=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_PE="$2"
      shift 2
      ;;
    # オプション-cut、--cut
    '-cut'|'--cut' )
      cut=$2
      # オプションに引数がなかった場合（必須）
      if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
        echo "$PROGNAME:「$1」オプションには引数が必要です" 1>&2
        exit 1
      fi
      ARG_CUT="$2"
      shift 2
      ;;
    '--'|'-' )
      # 「-」か「--」だけ打った時
      shift 1
      param+=( "$@" )
      break
      ;;
    -*)
      echo "$PROGNAME: 「$(echo $1 | sed 's/^-*//')」オプションは存在しません。'$PROGNAME -h'で確認してください" 1>&2
      exit 1
      ;;
    *)
      # コマンド引数（オプション以外のパラメータ）
      if [[ ! -z "$1" ]] && [[ ! "$1" =~ ^-+ ]]; then
        param+=( "$1" )
        shift 1
      fi
      ;;
  esac
done

# 「-d」のオプション指定がない場 
if [ -z $ARG_D ]; then
  echo "$PROGNAME:「-d」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-in」のオプション指定がない場 
if [ -z $ARG_IN ]; then
  echo "$PROGNAME:「-in」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-re」のオプション指定がない場 
if [ -z $ARG_RE ]; then
  echo "$PROGNAME:「-re」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-du」のオプション指定がない場 
if [ -z $ARG_DU ]; then
  echo "$PROGNAME:「-du」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-top」のオプション指定がない場 
if [ -z $ARG_TOP ]; then
  echo "$PROGNAME:「-top」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-it」のオプション指定がない場 
if [ -z $ARG_IT ]; then
  echo "$PROGNAME:「-it」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-cov」のオプション指定がない場 
if [ -z $ARG_COV ]; then
  echo "$PROGNAME:「-cov」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-int」のオプション指定がない場 
if [ -z $ARG_INT ]; then
  echo "$PROGNAME:「-int」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-pe」のオプション指定がない場 
if [ -z $ARG_PE ]; then
  echo "$PROGNAME:「-pe」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

# 「-cut」のオプション指定がない場 
if [ -z $ARG_CUT ]; then
  echo "$PROGNAME:「-cut」オプションは必須です。正しいオプションを指定してください" 1>&2
  echo $HELP_MSG 1>&2
  exit 1
fi

cd ${output_dir}
seqkit sort -lr ${reads_data} | seqkit head -n ${top} > top.fa
top_reads=`seqkit seq -ni top.fa`
/home/ozawat/local/src/minialign/minialign -x pacbio -f 0 -m 0.00001  top.fa ${reads_data}  -t16 | samtools view -Sb -F 4 | samtools sort > top_vs_reads_sorted.bam
touch repeat.fa
type > repeat.fa
for ref in ${top_reads}
do
  mkdir ${ref}
  cd ${ref}
  mkfifo hoge
  samtools view -h ../top_vs_reads_sorted.bam > hoge &
  STACK_YAML=${realigner_dir}/stack.yaml stack exec -- realigner -i ${it} ${ref} < hoge | samtools view -Sb | samtools sort >  ${ref}_vs_reads_realigned.bam
  samtools view -h ${ref}_vs_reads_realigned.bam | STACK_YAML=${dump_dir}/stack.yaml stack exec -- dump-consensus ${ref} -c > ${ref}_vs_reads_realigned_count.csv
  python ${src_dir}/count_to_consensus.py -s ${ref} -f ${ref}_vs_reads_realigned_count.csv -g True > consensus_${ref}_vs_reads_realigned.fa
  python ${src_dir}/bam-alignment_coverage.py ${ref}_vs_reads_realigned.bam > ${ref}_vs_reads_depth.txt
  python ${src_dir}/depth-high_coverage.py ${ref}_vs_reads_depth.txt  ${cov} > ${ref}_vs_reads_region.txt
  python ${src_dir}/bam-alignment_terminal.py ${ref}_vs_reads_realigned.bam> ${ref}_vs_reads_terminal.csv
  python ${src_dir}/terminal-find_peaks.py ${ref}_vs_reads_terminal.csv ${interval} ${pe} > ${ref}_vs_reads_peak.txt
  python ${src_dir}/peak-cut_region.py ${ref}_vs_reads_region.txt ${ref}_vs_reads_peak.txt ${cut} > ${ref}_vs_reads_cut.bed
  seqkit subseq --bed ${ref}_vs_reads_cut.bed consensus_${ref}_vs_reads_realigned.fa > ${ref}_vs_reads_result_gapped.fa
  python ${src_dir}/remove_gap_from_fasta.py ${ref}_vs_reads_result_gapped.fa> ${ref}_vs_reads_result.fa
  cd ..
  cat ./${ref}/${ref}_vs_reads_result.fa >> repeat.fa
done
python ${src_dir}/cal_cluster_size.py repeat.fa > cluster_size.csv
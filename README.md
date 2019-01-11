# make_rep_database

program of detecting and clustering repeat sequences from reads.

## Description

`make_rep_database` detects and classifies the repeat sequences from reads in FASTA file, and outputs the repeat databases with FASTA format.

## Requirement

- seqkit Version: 0.8.0
- samtools Version: 1.7 Using htslib 1.7
- minialign Version: 0.6.0-42-g4d6b274, Build: SSE4.1
- snakemake Version: 5.4.0

realignerとdump-consensusの設定は、ラボ内DocBase参照。

- realigner commit-number: 1b80e94e8139b6fa710876044651bc17b0405d06

- dump-consensus commit-number: a7733b5e0d89c7536f72c3f4e1d5a0c7747d013f

- multialn  comimt-number: 8102e5be6e1319fb07917c44cb484f564ddd91ff

- haskell-sam   commit-number: e8d1bc609eae5415755ca3084b308dcba12d4d7d

- fastlinkedlist    commit-number: bac809b90ed06b2fbb28b78da7c69021e3641d53

- bioseq    commit-number: d36162340f88169a137797124e437e44a6e688ac

- attoparsec-applicative    comimt-number: bd4b2902570dec6af0078171b94075f25f01e708

- python3

    - pysam　Version: 0.11.2.2

    - wiggelen Version: 0.4.1

    - scipy Version: 1.0.1

    - matplotlib Version: 2.2.2

    - biopython Version: 1.71

    - click Version: 7.0

## Usage

`make_rep_database.py -f input.fa -d output_directory -r realigner_directory -s dump-consensus_directory -t the_number_of_top_reads_to_use -i the_maximum_number_of_iterations_of_realigner -c threshold_of_depth -k threshold_of_peak_interval -p threshold_of_peak_coverage -l threshold_of_cut_interval`

今の仕様だと、

- `-f` 入力のFASTAファイル(絶対パス)

- `-d` 出力ファイルを出すディレクトリを指定(絶対パス)

- `-r` realignerのディレクトリ(絶対パス)

- `-s` dump-consensusのディレクトリ(絶対パス)

- `-t` 長い順に上から何本をrepeat検出のreferenceとして使うか

- `-i` realignerを最大何回回すか

- `-c` カバレッジいくつ以上の領域をリピートとみなすか

- `-k` 何bp内のリードの終端のカバレッジを極大値とみなすか

- `-p` リードの終端の数がいくつ以上のとき、ピークとみなすか

- `-l` リピート領域の切断位置をそれぞれ何bp以上離すか

(オプション多すぎるため、設定ファイルを用意したほうがよさそう。)

出力として、`-d`で指定したディレクトリに、topで使われた各リードの名前のフォルダ(中間ファイルが入っている), cluster_size.csv, top.fa, repeat.fa,top_vs_reads_sorted.bam, config.json, .snakemake(snakemakeのログ等がア￥入る)が入る。この中のrepeat.faがrepeatのデータベースとなる。今のところ、clusterのサイズは全て0.01になっている。

## Installation

    $ git clone https://github.com/heragaji/make_rep_database.git

## Author

Taro Ozawa

## License

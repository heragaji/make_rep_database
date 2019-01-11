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


## Usage

`snakemake --directory working_directory`

作業ディレクトリに、`config.json`という名前のJSONファイルを作成する。  
JSONファイルの中身は、

``` json
{
    "data": "入力のFASTAファイル(絶対パス)",
    "output": "出力ファイルを出すディレクトリを指定(絶対パス)",
    "realigner": "realignerのディレクトリ(絶対パス)",
    "dump_consensus": "dump-consensusのディレクトリ(絶対パス)",
    "top": "長い順に上から何本をrepeat検出のreferenceとして使うか",
    "iteration": "realignerを最大何回回すか",
    "coverage": "カバレッジいくつ以上の領域をリピートとみなすか",
    "interval": "何bp内のリードの終端のカバレッジを極大値とみなすか",
    "peak": "リードの終端の数がいくつ以上のとき、ピークとみなすか",
    "cut": "リピート領域の切断位置をそれぞれ何bp以上離すか",
    "src_dir": "Snakefileと同じ位置にある、srcフォルダを指定(絶対パス)"
}
```

出力として、"output"で指定したディレクトリに、topで使われた各リードの名前のフォルダ(中間ファイルが入っている), cluster_size.csv, top.fa, repeat.fa,top_vs_reads_sorted.bam, config.jsonが入る。この中のrepeat.faがrepeatのデータベースとなる。今のところ、clusterのサイズは全て0.01になっている。
作業ディレクトリには.snakemake(snakemakeのログ等が入る)が入る。

Snakefile自体のパスをプログラム内で得られれば、最後の"src_dir"はいらないので、なんとかしたい。
`--cores int`,`--qsub int`等でコア数やジョブの投入のオプションもできる。

## Installation

    $ git clone https://github.com/heragaji/make_rep_database.git

## Author

Taro Ozawa

## License

import os
import argparse
import shutil
import sys
import glob


def parse_arguments():
    parser = argparse.ArgumentParser(description="Pathseq report Script")
    parser.add_argument("--pathseq_dir", required=True, help="Path to the pathseq directory")
    parser.add_argument("--species", required=True, help="Species name")
    return parser.parse_args()


def copy_mapfile(pathseq_dir, dest_dir):
    mapfile_src = os.path.join(pathseq_dir, "mapfile")
    mapfile_dest = os.path.join(dest_dir, "mapfile")
    try:
        shutil.copy(mapfile_src, mapfile_dest)
        print(f"mapfile 已复制到 {mapfile_dest}")
    except Exception as e:
        print(f"复制 mapfile 时出错: {e}")
    return mapfile_dest


def create_run_sh(species):
    run_sh_content = f"""multi_pathseq \
 --mapfile mapfile \
 --genomeDir /SGRNJ06/randd/public/genome/rna/celescope_v2/{species} \
 --filter_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm10/mm10.fasta.img \
 --kmer_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_mm10/mm10.bfi \
 --microbe_bwa_image /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.fa.img \
 --microbe_dict /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_microbe.dict \
 --microbe_taxonomy_file /SGRNJ06/randd/USER/wangjingshen/pipeline/sc16S/reference/pathseq_microbe/pathseq_taxonomy.db \
 --thread 16 \
 --steps_run sample,starsolo,count_pathseq,analysis_pathseq"""
    with open("run.sh", "w") as f:
        f.write(run_sh_content)
    print("run.sh 已成功创建")


def process_mapfile(mapfile, pathseq_dir):
    with open(mapfile, "r") as f:
        lines = f.readlines()

    for line in lines:
        fields = line.strip().split()
        if len(fields) < 3:
            print(f"跳过格式不正确的行: {line}")
            continue

        sample = fields[2]
        sample_dir = os.path.join(os.getcwd(), f"{sample}/02.pathseq")
        os.makedirs(sample_dir, exist_ok=True)

        links = [
            (f"{pathseq_dir}/{sample}/02.pathseq/{sample}_pathseq.bam", f"{sample_dir}/{sample}_pathseq.bam"),
            (
                f"{pathseq_dir}/{sample}/02.pathseq/{sample}_pathseq_score.txt",
                f"{sample_dir}/{sample}_pathseq_score.txt",
            ),
            (f"{pathseq_dir}/{sample}/01.starsolo/outs/{sample}_*addRG.bam", f"{sample_dir}/{sample}_unmapped.bam"),
        ]

        for src, dst in links:
            src = glob.glob(src)[0] if "*" in src else src
            try:
                if not os.path.exists(src):
                    sys.exit(f"源文件不存在: {src}")
                if not os.path.exists(dst):
                    os.symlink(src, dst)
                    print(f"创建软链接: {dst}")
            except Exception as e:
                print(f"创建软链接时出错 ({src} -> {dst}): {e}")


def main():
    args = parse_arguments()
    pathseq_dir = args.pathseq_dir
    species = args.species

    # Step 1: 复制 mapfile
    mapfile = copy_mapfile(pathseq_dir, os.getcwd())

    # Step 2: 创建 run.sh
    create_run_sh(species)

    # Step 3: 处理 mapfile 并创建软链接
    process_mapfile(mapfile, pathseq_dir)


if __name__ == "__main__":
    main()

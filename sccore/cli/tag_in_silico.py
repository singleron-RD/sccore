import argparse
import random
from abc import abstractmethod
from collections import defaultdict

import pandas as pd
from scipy.stats import poisson


class Chip:
    def __init__(self, well_num, barcode_num, seed=1):
        self._well_num = well_num
        self._barcode_num = barcode_num
        self._seed = seed

        self._well_barcode = self.get_well_barcodes()

    def get_well_barcodes(self):
        """
        Set each well a random barcode
        Returns:
            well_barcode: {WELL index(0-based): BARCODE index(0-based)}
        """
        random.seed(self._seed)
        well_barcode = dict()
        for well in range(self._well_num):
            barcode = random.randint(0, self._barcode_num - 1)
            well_barcode[well] = barcode
        return well_barcode

    @property
    def well_num(self):
        return self._well_num

    @property
    def barcode_num(self):
        return self._barcode_num

    @property
    def well_barcode(self):
        return self._well_barcode


# model to estimate doublets
class Model:
    @abstractmethod
    def count(self, well_num: int, cell_num: int) -> tuple[int, int, int]:
        """
        Returns:
            number of wells with 0,1,2 cells; simplify the wells with >2 cells as 2 cells
        """
        pass


# Zero-Inflated Poisson
class ZIP_model(Model):
    def __init__(self, PSI: float = 0.4) -> None:
        self._PSI = PSI

    @staticmethod
    def ZIP_pmf(k: int, mu: float, PSI: float) -> float:
        """
        Zero-Inflated Poisson Probability Mass Function

        Args:
            k: The number of occurrences.
            mu: The mean of the Poisson distribution.
            PSI: The zero-inflation parameter (probability of extra zeros).
        Returns:
            The probability of observing exactly `k` occurrences under the ZIP distribution.
        """
        if k == 0:
            return float(PSI + (1 - PSI) * poisson.pmf(0, mu=mu))
        else:
            return float((1 - PSI) * poisson.pmf(k, mu=mu))

    def count(self, well_num: int, cell_num: int) -> tuple[int, int, int]:
        mu = cell_num / well_num
        count0 = int(self.ZIP_pmf(k=0, mu=mu, PSI=self._PSI) * well_num)
        count1 = int(self.ZIP_pmf(k=1, mu=mu, PSI=self._PSI) * well_num)
        count2 = well_num - count0 - count1
        return count0, count1, count2


# Simulation, assume cells are uniformly distributed across tags
class AddCell:
    def __init__(
        self,
        cell_num: int,
        tag_num: int,
        chip: Chip,
        model: Model,
        seed=1,
    ):
        self._cell_num = cell_num
        self._tag_num = tag_num
        self._chip = chip
        self._model = model
        self._seed = seed

    def _get_well_cellNum(self) -> dict[int, int]:
        """
        Returns:
            well_cellNum: {well_index: cell number}
        """
        random.seed(self._seed)
        well_cellNum = dict()
        well_num = self._chip.well_num
        count0, count1, count2 = self._model.count(well_num, self._cell_num)
        wells = list(range(well_num))
        random.shuffle(wells)

        for w in wells[:count0]:
            well_cellNum[w] = 0
        for w in wells[count0 : count0 + count1]:
            well_cellNum[w] = 1
        for w in wells[count0 + count1 : well_num]:
            well_cellNum[w] = 2
        return well_cellNum

    def _get_well_tag(self, well_cellNum) -> dict[int, list]:
        """
        Returns:
            tag_dict: {well_index: tag}
        """
        random.seed(self._seed)
        tag_dict = {}
        for well in well_cellNum:
            if well_cellNum[well] == 2:
                tag_dict[well] = [random.randint(0, self._tag_num - 1), random.randint(0, self._tag_num - 1)]
            elif well_cellNum[well] == 1:
                tag_dict[well] = [random.randint(0, self._tag_num - 1)]
        return tag_dict

    def summarize(self, well_tag):
        single = 0
        classified = 0
        class_dict = defaultdict(lambda: defaultdict(int))
        for well in well_tag:
            barcode = self._chip.well_barcode[well]
            for tag in well_tag[well]:
                class_dict[barcode][tag] += 1
        for barcode in class_dict:
            tag_number = sum(class_dict[barcode].values())
            if tag_number > 1:
                class_flag = True
                for tag in class_dict[barcode]:
                    if class_dict[barcode][tag] > 1:
                        class_flag = False
                        break
                if class_flag:
                    classified += 1
            else:
                single += 1
        recovered = len(class_dict)
        unclassified = recovered - single - classified
        # print(recovered, single, classified, unclassified)
        return recovered, single, classified, unclassified

    def run_once(self):
        well_cellNum = self._get_well_cellNum()
        well_tag = self._get_well_tag(well_cellNum)
        return self.summarize(well_tag)

    def run(self, iter: int) -> tuple[int, int, int, int]:
        """
        Mean of the iteration results.
        Args:
            iter: int, number of iterations.
        Returns:
            tuple: (recovered, single, classified, unclassified)
                - recovered: recovered cells.
                - single: the number of cells with tags == 1.
                - classified: the number of cells which have tags > 1 and the tags are not the same.
                - unclassified: the number of cells which have tags > 1 and the tags are the same.
        """
        results = [self.run_once() for _ in range(iter)]
        print(results)
        recovered, single, classified, unclassified = [int(sum(x) / iter) for x in zip(*results)]
        return recovered, single, classified, unclassified

    @property
    def cell_num(self):
        """number of loaded cells"""
        return self._cell_num

    @property
    def tag_num(self):
        return self._tag_num


def format_percent(number):
    return round(number * 100, 2)


def main():
    parser = argparse.ArgumentParser(description="clindex simulation")
    parser.add_argument("--cell", required=True, type=str)
    parser.add_argument("--well", required=True, type=int)
    parser.add_argument("--tag", required=True, type=str)
    parser.add_argument("--barcode_num", default=96 * 96 * 192, type=int)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    chip = Chip(well_num=args.well, barcode_num=args.barcode_num, seed=args.seed)
    model = ZIP_model(PSI=0.4)

    rows = []
    for cell_num in args.cell.split(","):
        cell_num = int(cell_num)
        for tag_num in args.tag.split(","):
            tag_num = int(tag_num)
            add_cell = AddCell(
                cell_num=cell_num,
                tag_num=tag_num,
                chip=chip,
                model=model,
                seed=args.seed,
            )
            recovered, single, classified, unclassified = add_cell.run(iter=2)
            df_dict = {
                "Cells loaded": cell_num,
                "Estimated recovered cells": recovered,
                "N_tags": tag_num,
                "Estimated doublet percent": format_percent(1 - single / recovered),
                "Estimated classified percent": format_percent(classified / recovered),
                "Estimated unclassified percent": format_percent(unclassified / recovered),
            }
            rows.append(df_dict)
    df = pd.DataFrame(rows)
    print(df)
    df.to_csv(f"well{args.well}_tag{args.tag}_seed{args.seed}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

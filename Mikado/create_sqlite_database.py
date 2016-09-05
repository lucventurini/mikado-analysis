#!/usr/bin/env python3

import sqlalchemy
from sqlalchemy import Column, ForeignKey, PrimaryKeyConstraint
from sqlalchemy.types import Integer, LargeBinary, String, Boolean, Float
from sqlalchemy.engine import create_engine
from sqlalchemy.orm.session import sessionmaker
import sqlite3
from Mikado.utilities.log_utils import create_default_logger
import argparse
from sqlalchemy.engine.reflection import Inspector
from sqlalchemy.ext.declarative import declarative_base
import os
import functools



DBBASE = declarative_base()
__doc__ = """"Quick script to load into a database the comparison data for Mikado, as BLOBs"""


class Indexer(DBBASE):

    __tablename__ = "indexer"
    __table_args__ = {"extend_existing": True}

    m_index = Column(Integer, primary_key=True)
    species = Column(String(100), unique=False)
    aligner = Column(String(100), unique=False)
    assembler = Column(String(100), unique=False)
    constraint = PrimaryKeyConstraint("species", "aligner", "method",
                                      name="method_combination")

    def __init__(self, species, aligner, assembler):

        self.species, self.aligner, self.assembler = species, aligner, assembler


class CompareFiles(DBBASE):

    __tablename__ = "compare_files"
    __maxsize__ = 40*1024*1024

    m_index = Column(Integer, ForeignKey(Indexer.m_index), unique=False, nullable=False)
    filtered = Column(Boolean, unique=False)
    constraint = PrimaryKeyConstraint("m_index", "filtered", name="compare_idx")
    refmap = Column(LargeBinary(length=__maxsize__))
    tmap = Column(LargeBinary(length=__maxsize__))
    stats = Column(LargeBinary(length=__maxsize__))
    __table_args__ = (constraint,)

    def __init__(self, index, base, filtered):

        if not isinstance(filtered, bool):
            raise ValueError("Invalid filtered value: {}")

        assert index is not None
        self.m_index = index
        self.filtered = filtered
        self.tmap = memoryview(open("{}.tmap".format(base), "rb").read())
        self.refmap = memoryview(open("{}.refmap".format(base), "rb").read())
        self.stats = memoryview(open("{}.stats".format(base), "rb").read())
        print(self.m_index, base, self.filtered)


class Statistics(DBBASE):

    __tablename__ = "basic_statistics"
    __table_args__ = {"extend_existing": True}

    id = Column(Integer, ForeignKey(Indexer.m_index), primary_key=True)
    level = Column(String(20))
    constraint = PrimaryKeyConstraint("id", "level", name="stat_combination")
    precision = Column(Float)
    recall = Column(Float)
    f1 = Column(Float)


def main():

    logger = create_default_logger("stat_serializer")

    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--db", required=True, help="SQLite database to connect to.")
    parser.add_argument("--force", action="store_true", default=False)
    parser.add_argument("input_files",
                        help="""TXT tab-delimited file, specifying the input files in the following way:
                        - species
                        - aligner
                        - assembler
                        - basename for the comparisons against the complete reference
                        - basename for the comparisons against the filtered reference
                        """)
    args = parser.parse_args()


    # Create the database
    connector = functools.partial(sqlite3.connect,
                                  database=args.db, check_same_thread=False)
    engine = create_engine("sqlite://", creator=connector)

    if args.force is True:
        logger.warn("Removing old data because force option in place")
        meta = sqlalchemy.MetaData(bind=engine)
        meta.reflect(engine)
        for tab in reversed(meta.sorted_tables):
            logger.warn("Dropping %s", tab)
            tab.drop()

    inspector = Inspector.from_engine(engine)
    Session = sessionmaker(bind=engine)
    # session = Session(bind=engine, autocommit=True, autoflush=True)
    session = Session()

    if Indexer.__tablename__ not in inspector.get_table_names():
        DBBASE.metadata.create_all(engine)  # @UndefinedVariable


    with open(args.input_files) as input_files:

        for row in input_files:
            species, aligner, assembler, complete, filtered = row.rstrip().split()
            if not os.path.exists("{}.stats".format(complete)):
                raise ValueError("Original file not found; line:\n{}".format(row))
            if not os.path.exists("{}.stats".format(filtered)):
                raise ValueError("Filtered file {} not found; line:\n{}".format(
                    "{}.stats".format(filtered),
                    row))
            current_species = Indexer(species, aligner, assembler)
            session.add(current_species)
            session.commit()
            print(current_species.m_index)
            complete_load = CompareFiles(current_species.m_index, complete, filtered=False)
            session.add(complete_load)
            # session.commit()
            # filtered_load = CompareFiles(current_species.index, filtered, filtered=True)
            # session.add(filtered_load)
            # session.commit()
            #
            # orig_lines = [line.rstrip() for line in open(orig)]
            # filtered_lines = [line.rstrip() for line in open(filtered)]
            # # In the stats we have precision as second and sensitivity as first,
            # # we have to invert
            # for index, line_index in enumerate([5, 7, 8, 9, 12, 15]):
            #     precision = float(orig_lines[line_index].split(":")[1].split()[1])
            #     recall = float(filtered_lines[line_index].split(":")[1].split()[0])
            #     stats[
            #         # Name of the statistic:Base, Exon, etc
            #         list(stats.keys())[index]][
            #         b"TopHat"].append((precision, recall))
    session.commit()

main()
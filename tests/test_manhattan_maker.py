import pytest
import argparse
import sys
import os
import pandas as pd

import matplotlib.pyplot as plt

from manhattan_maker import parse_args, check_args, read_input_file, create_manhattan_plot
from manhattan_maker import ProgramError


# Testing the manhattan_maker.py file

########################################################################################################################
# Testing the parse_args function

def test_parse_args():
    """Test the parse_args function. Make sure it returns an argparse.Namespace object"""
    print('CONTENTS OF SYS.ARGV: ', sys.argv)
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(args=['--in-file', filename])

    assert type(args) == argparse.Namespace


def test_parse_args_filename():
    """Test the parse_args function. Make sure the filename is correct"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(args=['--in-file', filename])

    assert args.in_file == filename


def test_parse_args_column_names():
    """Test the parse_args function. Make sure the column names change if the user specifies them"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(
        args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos', 'position',
              '--col-pvalue', 'pval'])

    assert args.col_name == 'rsid'
    assert args.col_chr == 'chromosome'
    assert args.col_pos == 'position'
    assert args.col_pvalue == 'pval'


def test_parse_args_exclude_chr():
    """Test the parse_args function. Make sure the chromosome is excluded if the user specifies it"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(
        args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos', 'position',
              '--col-pvalue', 'pval', '--exclude-chr', '1'])

    assert args.exclude_chr == '1'


def test_parse_args_only_chr():
    """Test the parse_args function. Make sure the chromosome is excluded if the user specifies it"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(
        args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos', 'position',
              '--col-pvalue', 'pval', '--only-chr', '1'])

    assert args.only_chr == '1'


def test_parse_args_groups():
    """Test the parse_args function. Make sure the groups are correctly read"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(
        args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos', 'position',
              '--col-pvalue', 'pval', '--use-groups', '--group-colors', 'group1,group2'])

    assert args.use_groups_flag == True
    assert args.group_colors == 'group1,group2'


########################################################################################################################
# Testing the check_args function

def test_check_args_lacking_input_file():
    """Test the check_args function. Make sure 'check_args' rises an exception if the input file is not specified"""
    args = parse_args(args=['--output', 'test_output.txt'])
    with pytest.raises(Exception) as e_info:  # e_info saves the exception object, so you can extract details from it.
        check_args(args)


def test_check_args_exclude_chr_and_only_chr():
    """Test the check_args function. Make sure 'check_args' raises an exception if the user specifies both
    exclude_chr and only_chr """
    args = parse_args(
        args=['--in-file', 'GWAS_output_testfile.txt', '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
              'position', '--col-pvalue', 'pval', '--exclude-chr', '1', '--only-chr', '2'])
    with pytest.raises(
            ProgramError) as e_info:  # e_info saves the exception object, so you can extract details from it.
        check_args(args)


def test_check_args_group_colors_without_use_groups():
    """Test the check_args function. Make sure 'check_args' raises an exception if the user specifies group_colors
    without use_groups"""
    args = parse_args(
        args=['--in-file', 'GWAS_output_testfile.txt', '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
              'position', '--col-pvalue', 'pval', '--group-colors', 'group1,group2,group3'])
    # TODO: Deal with two color groups case
    with pytest.raises(ProgramError) as e_info:  # e_info saves the exception object, so you can extract details
        check_args(args)


########################################################################################################################
# Testing the read_input_file function

def test_read_input_file():
    """Test the read_input_file function. Make sure it reads the file fine"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval'])
    df = read_input_file(filename, args.use_groups_flag, args)
    print(df)
    assert type(df) == pd.DataFrame
    assert df.columns.tolist() == ['chrom', 'pos', 'snp', 'conf']


def test_read_input_file_missing_columns():
    """Test the read_input_file function. Make sure it reads the file fine"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position'])
    with pytest.raises(ProgramError) as e_info:
        read_input_file(filename, True, args)


def test_read_input_file_use_groups():
    """Test the read_input_file function. Make sure it reads the file fine"""
    filename = 'GWAS_output_testfile.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--use-groups', '--col-group', 'group'])
    check_args(args)
    df = read_input_file(filename, args.use_groups_flag, args)
    print(df)
    assert type(df) == pd.DataFrame
    assert df.columns.tolist() == ['chrom', 'pos', 'snp', 'conf', 'group']


def test_read_input_file_no_group_column():
    """Test the read_input_file function. Make sure it reads the file fine"""
    filename = 'GWAS_output_testfile_2.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--use-groups'])
    with pytest.raises(ProgramError) as e_info:
        read_input_file(filename, args.use_groups_flag, args)


def test_read_input_file_incomplete_group_column():
    """Test the read_input_file function. Make sure it reads the file fine even if the 'group' column is incomplete"""
    filename = 'GWAS_output_testfile_3.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--col-group', 'group', '--use-groups',
                            '--group-colors', '#ff0000,#00ff00,#0000ff'])
    df = read_input_file(filename, args.use_groups_flag, args)
    assert type(df) == pd.DataFrame
    assert df.columns.tolist() == ['chrom', 'pos', 'snp', 'conf', 'group']


def test_read_input_file_not_enough_colors():
    """Test the read_input_file function. Make sure it reads the file fine"""
    filename = 'GWAS_output_testfile_3.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--col-group', 'group', '--use-groups',
                            '--group-colors', '#ff0000,#00ff00'])
    with pytest.raises(ProgramError) as e_info:
        read_input_file(filename, args.use_groups_flag, args)

    #TODO: evidence of this test passing when it shouldn't be.


########################################################################################################################
# Testing the create_manhattan_plot function

def test_create_manhattan_plot(monkeypatch):
    """Test the create_manhattan_plot function. Make sure it creates the manhattan plot fine"""
    filename = 'GWAS_output_testfile.txt'
    # Remove file 'manhattan.png' if it exists
    try:
        os.remove('manhattan.png')
    except OSError:
        pass
    # Confirm there is no file named 'manhattan.png'
    assert not os.path.exists('manhattan.png')

    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval'])
    df = read_input_file(filename, args.use_groups_flag, args)
    monkeypatch.setattr(plt, 'show', lambda: None)  # This is to avoid the plot being shown on the screen
    create_manhattan_plot(df, args)
    assert os.path.isfile('manhattan.png')  # Check if the file has been created


def test_create_manhattan_plot_with_groups(monkeypatch):
    """Test the create_manhattan_plot function. Make sure it creates the manhattan plot fine"""
    filename = 'GWAS_output_testfile.txt'
    # Remove file 'manhattan.png' if it exists
    try:
        os.remove('manhattan.png')
    except OSError:
        pass
    # Confirm there is no file named 'manhattan.png'
    assert not os.path.exists('manhattan.png')

    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--use-groups', '--col-group', 'group'])
    check_args(args)
    df = read_input_file(filename, args.use_groups_flag, args)
    monkeypatch.setattr(plt, 'show', lambda: None)  # This is to avoid the plot being shown on the screen
    create_manhattan_plot(df, args)
    assert os.path.isfile('manhattan.png')  # Check if the file has been created


def test_create_manhattan_plot_given_colors(monkeypatch):
    """Test the create_manhattan_plot function. Make sure it creates the manhattan plot fine"""
    filename = 'GWAS_output_testfile_3.txt'
    # Remove file 'manhattan.png' if it exists
    try:
        os.remove('manhattan.png')
    except OSError:
        pass
    # Confirm there is no file named 'manhattan.png'
    assert not os.path.exists('manhattan.png')

    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--significant-color', '#808080', '--use-groups',
                            '--col-group', 'group', '--group-colors', '#ff0000,#00ff00,#0000ff'])
    color_list = args.group_colors.split(',')
    assert len(color_list) == 3

    df = read_input_file(filename, args.use_groups_flag, args)
    monkeypatch.setattr(plt, 'show', lambda: None)  # This is to avoid the plot being shown on the screen
    create_manhattan_plot(df, args)
    assert os.path.isfile('manhattan.png')  # Check if the file has been created


########################################################################################################################
# Integration tests

def test_parser_parsercheck_reader():

    filename = 'GWAS_output_testfile_3.txt'
    args = parse_args(args=['--in-file', filename, '--col-name', 'rsid', '--col-chr', 'chromosome', '--col-pos',
                            'position', '--col-pvalue', 'pval', '--significant-color', '#808080', '--use-groups',
                            '--col-group', 'group', '--group-colors', '#ff0000,#00ff00,#0000ff'])
    check_args(args)
    df = read_input_file(filename, args.use_groups_flag, args)
    assert type(df) == pd.DataFrame
    assert df.columns.tolist() == ['chrom', 'pos', 'snp', 'conf', 'group']
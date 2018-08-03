##
# \file multicomparsion_parse
# \brief Class for handling Multicomparison's "all_robustness" output files.
#

from copy import deepcopy


class DataParser:

    def __init__(self):
        ## Matrix for the log.
        self.matrix = []
        self.headers = {}

    ##
    # \fn copy
    # \brief Perform a deep copy of the DataParser.
    # @return Returns a
    def copy(self):
        cop = DataParser()
        cop.matrix = self.copy_matrix()
        cop.headers = self.copy_headers()
        return cop

    ##
    # \fn copy_headers
    # \brief Copy the header dictionary of this DataParser.
    # @return Returns the copied header dictionary.
    def copy_headers(self):
        return deepcopy(self.headers)

    ##
    # \fn copy_matrix
    # \brief Copy the matrix of this DataParser
    # @return Returns a deep copy of the internal matrix.
    def copy_matrix(self):
        return deepcopy(self.matrix)

    ##
    # \fn copy_headless_matrix
    # \brief Copies the internal matrix and returns it headless.
    # @return Returns a deep copy of the internal matrix, but without the first
    #     row.
    def copy_headless_matrix(self):
        return deepcopy(self.matrix[1:])

    ##
    # \fn copy_headers_only
    # \brief Copy the data parser, but with no data stored.
    def copy_headers_only(self):
        cop = DataParser()
        cop.headers = self.copy_headers()
        cop.matrix = [deepcopy(self.matrix[0])]
        return cop

    ##
    # \fn rename_headers
    # \brief Renames the headers and returns it.
    # @param rename_dict A dictionary that has keys be the headers of the
    #   current matrix is, and the values be a string of the new headers.
    def rename_headers(self, rename_dict):
        cop = self.copy()
        cop.headers = {}
        new_first_row = cop.matrix[0]
        try:
            for key in self.headers.keys():
                # Set the headers.
                cop.headers[rename_dict[key]] = self.headers[key]
                # Set the first row's values.
                cop.matrix[0][self.headers[key]] = rename_dict[key]
            return cop
        except KeyError:
            raise KeyError("rename_dict was missing header key "+key)

    ##
    # \fn load_matrix
    # \brief Load the Multicomparsion output into this data parser.
    #
    # @param fp Filepath to read in.
    # @param delim The delimiter for the input data. By default, it will match
    #         all white space.
    def load_matrix(self, fp, delim=None):
        self.matrix = []
        with open(fp, 'r') as f:
            lines = f.readlines()
            # Append headers first
            self.matrix.append(lines[0].strip().split(delim))
            # Iterate through the rest of the data.
            for col_str in lines[1:]:
                col_list = col_str.strip().split(delim)
                self.matrix.append(list(map(float, col_list)))

        for i, col in enumerate(self.matrix[0]):
            self.headers[col] = i

    ##
    # \fn to_csv
    # \brief Convert the weird format of the multiComparison output to CSV.
    # @param outpath Print to a file?
    # @param include_headers Output with headers.
    # @return If no outpath specified, will return the output as a string.
    #     Else returns None
    def to_csv(self, outpath=None, include_headers=True):
        output = ""

        starting_index = 0
        if not include_headers:
            starting_index = 1

        for row in self.matrix[starting_index:]:
            output += ",".join(list(map(str, row))) + "\n"
        if outpath is None:
            # Slice off the newline at the end.
            return output[:-1]
        else:
            with open(outpath, 'w') as f:
                f.write(output)
            return None

    ##
    # \fn get_zero_filtered
    # \brief Returns the line entries that do not contain all zeros for each of
    #     the given items in the filter_list
    # @param filter_list If the entry has only zeros in the columns indicated,
    #     remove them from the DataParser.
    def get_zero_filtered(
            self, filter_list=["early", "stat", "n_stat", "dyn"]):
        new_parser = DataParser()
        # Begin the matrix with nothing but the headers.
        new_parser.headers = self.copy_headers()
        new_matrix = [deepcopy(self.matrix[0])]
        for row in self.matrix[1:]:
            use_row = False
            for head in filter_list:
                if head in self.headers:
                    if row[self.headers[head]] != 0:
                        use_row = True
            if use_row:
                new_matrix.append(deepcopy(row))

        new_parser.matrix = new_matrix

        return new_parser

    ##
    # \fn get_filtered
    # \brief Returns a new parser object, but only has entries that have the
    #     given key equal to the specified value.
    def get_filtered(self, key, value):
        new_parser = DataParser()
        # Begin the matrix with nothing but the headers.
        new_parser.headers = self.copy_headers()
        new_matrix = [deepcopy(self.matrix[0])]
        for row in self.matrix[1:]:
            use_row = False
            if key in self.headers:
                    if row[self.headers[key]] == value:
                        use_row = True
            if use_row:
                new_matrix.append(deepcopy(row))

        new_parser.matrix = new_matrix
        return new_parser

    ##
    # \fn get_filtered_less
    # \brief Returns a new parser object, but only has entries that have the
    #     given key less or equal to than another key's value.
    def get_filtered_less(self, less_key, greater_key):
        new_parser = DataParser()
        # Begin the matrix with nothing but the headers.
        new_parser.headers = self.copy_headers()
        new_matrix = [deepcopy(self.matrix[0])]
        for row in self.matrix[1:]:
            use_row = False
            if row[self.headers[less_key]] <= row[self.headers[greater_key]]:
                    use_row = True
            if use_row:
                new_matrix.append(deepcopy(row))
        new_parser.matrix = new_matrix
        return new_parser

    ##
    # \fn get_means
    # \brief Retrieve a new parser that only has a header row and an average
    #     row.
    def get_means(self):
        average_parser = DataParser()
        average_parser.headers = self.copy_headers()

        sum_row = [0]*len(self.matrix[0])
        for row in self.matrix[1:]:
            sum_row = [sum_row[i] + row[i] for i in range(len(row))]
        average_row = [float(e)/(len(self.matrix)-1) for e in sum_row]

        average_parser.matrix = [deepcopy(self.matrix[0]), average_row]
        return average_parser

    ##
    # \fn get_groupings
    # \brief Retrieve a List of DataParsers that each contain exactly one
    #     grouping from the given grouping keys.
    #
    # @param keys List of strings indicating which columns we should group by.
    # @return Returns a dict of the form {group_tuple : DataParser}.
    def get_groupings(self, keys):
        if type(keys) != list:
            raise TypeError("Was expecting keys to be a list.")

        # Dictionary to store all the groupings.
        groupings = {}
        for row in self.matrix[1:]:
            group_tuple = ()  # Yes, this is an empty tuple.

            for k in keys:
                # Make the group tuple from the values in the
                group_tuple += (row[self.headers[k]],)
            # We now have the correct tuple to check for this row

            # Check whether the grouping has been made yet.
            if group_tuple in groupings:
                # This grouping already exists, so add this row to that
                # grouping.
                groupings[group_tuple].matrix.append(deepcopy(row))
            else:
                # This grouping didn't yet exist, so create it.
                groupings[group_tuple] = DataParser()
                groupings[group_tuple].headers = self.copy_headers()
                groupings[group_tuple].matrix = (
                    [deepcopy(self.matrix[0])] + [deepcopy(row)])
        return groupings

    ##
    # \fn append_entries
    # \brief Append the contents of another parser to this parser.
    # \note Note that the otherparser must have the same headers as this
    #       parser.
    #
    # @param otherparser The parser to retrieve entries from.
    def append_entries(self, otherparser):
        append_parser = DataParser()
        append_parser.headers = self.copy_headers()
        append_parser.matrix = self.copy_matrix()

        for h in self.matrix[0]:
            if not h in otherparser.headers:
                raise RuntimeError(
                    "Could not append entries; headers did not match")

        # TODO:
        # Make this smarter. It's blindly appending. Instead, smart append
        # so that the order of the columns do not matter.
        for row in otherparser.copy_headless_matrix():
            append_parser.matrix.append(row)

        return append_parser

    ##
    # \fn get_stddevs
    # \brief Retrieve a parser that has each of the stand
    def get_stddevs(self):
        average_parser = self.get_means()
        stddev_parser = DataParser()
        stddev_parser.headers = self.copy_headers()
        stddev_parser.matrix = [deepcopy(self.matrix[0])]

        # Make a matrix formed of the difference to the mean squared.
        sqdiff_matrix = []
        for row in self.matrix[1:]:
            sqdiff_row = []
            for i, dat in enumerate(row):
                diff = dat - average_parser.matrix[1][i]
                # Calculate the sum squared
                sqdiff = (diff**2)
                sqdiff_row.append(sqdiff)
            sqdiff_matrix.append(sqdiff_row)

        # Get a summation row
        sum_row = [0]*len(sqdiff_matrix[0])
        for row in sqdiff_matrix:
            for i in range(len(row)):
                sum_row[i] = sum_row[i] + row[i]

        # Finally, calculate the standard deviation row
        stddev_row = [
            (col/(self.get_entry_count()-1))**(0.5) for col in sum_row]
        stddev_parser.matrix.append(stddev_row)

        return stddev_parser

    ##
    # \fn get_col_stderr
    # \brief Retrieve a column's standard error.
    def get_col_stderr(self, key):
        stddev_parser = self.get_stddevs()
        return (
            stddev_parser.matrix[1][stddev_parser.headers[key]]
            / (self.get_entry_count()**(0.5))
        )

    ##
    # \fn get_90_confint
    # \brief Retrieve a column's 90% confidence interval for the mean.
    # \note TODO: This really needs to be fixed, but only 2 days left for
    #     research. Hard coding in the value now.
    def get_90_confint(self, key):
        stderr = self.get_col_stderr(key)
        return stderr*1.645

    # -------------------------------------------------------------------------

    ##
    # \fn get_entry_count
    # \brief Return the number of data rows in the matrix.
    def get_entry_count(self):
        return len(self.matrix) - 1

    ##
    # \fn sort_on_key
    # \brief Sort the data based on a key
    # @param key Key (column label) to sort with.
    def sort_on_key(self, key):
        sorted_parser = DataParser()
        sorted_parser.headers = deepcopy(self.headers)
        sorted_parser.matrix = [deepcopy(self.matrix[0])]

        to_sort = self.matrix[1:]
        is_sorted = list(sorted(to_sort, key=(lambda r: r[self.headers[key]])))
        sorted_parser.matrix += is_sorted
        return sorted_parser

    ##
    # \fn only_columns
    # \brief Retrieve a new parser with only the specified keys, in order.
    #
    # @param keys Column keys to use.
    # @return Returns a new parser object with only the specified columns.
    def only_columns(self, keys):
        new_parser = DataParser()
        # Start the matrix off with the given keys
        new_parser.matrix = [deepcopy(keys)]
        # Rewrite the headers
        for i, k in enumerate(keys):
            new_parser.headers[k] = i

        for row in self.matrix[1:]:
            new_entry = [0]*len(keys)
            for k in keys:
                # Pick the element in the row with the right key,
                # and assign it to the new element in the new_entry
                new_entry[new_parser.headers[k]] = row[self.headers[k]]
            new_parser.matrix.append(new_entry)
        return new_parser

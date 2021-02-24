from caravagene import ConstructList
my_constructs = ConstructList("crecombio_spreadsheet.xlsx",
                              font='Ubuntu',
                              width=1200)
# my_constructs.to_pdf('crecombio_plot.pdf')
my_constructs.to_image('crecombio_plot.jpg')
from flask import Flask, request, render_template
from flask_wtf import FlaskForm
from wtforms.fields import StringField, SelectField
from wtforms.validators import DataRequired, ValidationError, regexp, Length, NumberRange
from lifelines import KaplanMeierFitter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pdb


app = Flask(__name__)
app.config['CSPR_ENABLE'] = True
app.config['SECRET_KEY'] = 'you-will-never-guess'

XdataSetX = [('GSE12417', 'GSE12417'), ('GSE71014',
                                        'GSE71014'), ('GSE22762', 'GSE22762')]


class SingleForm(FlaskForm):

    DataBase = SelectField("DataBase", choices=XdataSetX)
    GeneName = StringField(
        'GeneName', validators=[
            DataRequired(
                message='there must be input')])
    Low = StringField('Low', validators=[DataRequired()])
    High = StringField('High', validators=[DataRequired()])


class DoubleForm(FlaskForm):
    DataBase = SelectField("DataBase", choices=XdataSetX)
    GeneName1 = StringField('GeneName1', validators=[DataRequired()])
    GeneName2 = StringField('GeneName2', validators=[DataRequired()])
    Low1 = StringField('Low1', validators=[DataRequired()])
    High1 = StringField('High1', validators=[DataRequired()])
    Low2 = StringField('Low2', validators=[DataRequired()])
    High2 = StringField('High2', validators=[DataRequired()])


def ReadData(database, Gene):
    with open('static/' + database + '.csv', 'r') as x:
        os = np.recfromtxt('static/' + database + ".os")
        exp = [line.split(",")[2:]
               for line in x if Gene in line.split(",")[1].split(" /// ")]
        exp = [np.array(x, dtype='float32') for x in exp]  # turn to numpy
        exp = sum(exp) / len(exp)  # should it be average value?
        # expression, day, status
        data = sorted(zip(exp, os[:, 0], os[:, 1], range(0, len(os) + 1)))
        data = np.array((data))
        return data, os, exp.mean(), exp.std()


def single_submit(form):
    if form.validate_on_submit():

        database = form.DataBase.data
        Gene = form.GeneName.data
        low = int(form.Low.data)
        high = int(form.High.data)

        static = {}
        data, os, static['mean'], static['std'] = ReadData(database, Gene)

        num = len(os)
        low = max(int(num * low / 100), 1)
        high = max(int(num * high / 100), 1)

        Low, High = data[:, 1][0:low], data[:, 1][-high:]
        group1, group2 = data[:, 2][0:low], data[:, 2][-high:]

        kmf = KaplanMeierFitter()
        kmf.fit(Low, group1, label=Gene + '/low')
        ax = kmf.plot()
        kmf.fit(High, group2, label=Gene + '/high')
        kmf.plot(ax=ax)
        plt.savefig("static/test.png", bbox_inches='tight')
        plt.close()

        return render_template("single.html",
                               form=form,
                               image="test.png",
                               refresh=np.random.randn(),
                               static=static)
    else:
        return render_template("single.html",
                               form=form,
                               err=form.errors)


def double_submit(form):
    if form.validate_on_submit():

        database = form.DataBase.data
        Gene1, low1, high1 = form.GeneName1.data, int(
            form.Low1.data), int(form.High1.data)
        Gene2, low2, high2 = form.GeneName2.data, int(
            form.Low2.data), int(form.High2.data)

        static = {}
        data1, os, static['mean'], static['std'] = ReadData(database, Gene1)
        data2, os, static['mean'], static['std'] = ReadData(database, Gene2)

        num = len(os)
        low1, high1 = max(int(num * low1 / 100),
                          1), max(int(num * high1 / 100), 1)
        low2, high2 = max(int(num * low2 / 100),
                          1), max(int(num * high2 / 100), 1)

        low_low = np.array(list(set(data1[:, 3][0:low1]) & set(
            data2[:, 3][0:low2])), dtype='int32')
        low_high = np.array(list(set(data1[:, 3][0:low1]) & set(
            data2[:, 3][-high2:])), dtype='int32')
        high_low = np.array(
            list(set(data1[:, 3][-high1:]) & set(data2[:, 3][0:low2])), dtype='int32')
        high_high = np.array(
            list(set(data1[:, 3][-high1:]) & set(data2[:, 3][-high2:])), dtype='int32')

        # os day stauts

        kmf = KaplanMeierFitter()
        kmf.fit(os[:, 0][low_low], os[:, 1][low_low],
                label=Gene1 + '/low-' + Gene2 + '/low')
        ax = kmf.plot()
        kmf.fit(os[:, 0][low_high], os[:, 1][low_high],
                label=Gene1 + '/low-' + Gene2 + '/high')
        kmf.plot(ax=ax)
        kmf.fit(os[:, 0][high_low], os[:, 1][high_low],
                label=Gene1 + '/high-' + Gene2 + '/low')
        kmf.plot(ax=ax)
        kmf.fit(os[:, 0][high_high], os[:, 1][high_high],
                label=Gene1 + '/high-' + Gene2 + '/high')
        kmf.plot(ax=ax)

        plt.savefig("static/double.png", bbox_inches='tight')
        plt.close()

        return render_template("double.html",
                               form=form,
                               image="double.png",
                               refresh=np.random.randn(),
                               static=static)
    else:
        return render_template("double.html",
                               form=form)


@app.route("/single", methods=['GET', 'POST'])
def single_page():

    form = SingleForm(request.form)
    if form.is_submitted():
        return single_submit(form)
    return render_template("single.html",
                           form=form)


@app.route("/double", methods=['GET', 'POST'])
def double_page():

    form = DoubleForm(request.form)
    if form.is_submitted():
        return double_submit(form)
    return render_template("double.html", form=form)


@app.route("/index.html")
@app.route("/")
def index():
    return render_template("index.html")


@app.route("/test.html")
def test():
    return render_template("test.html")


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')

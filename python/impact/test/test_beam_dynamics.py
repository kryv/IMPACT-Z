import numpy as np
import unittest

from numpy.testing import assert_almost_equal as assert_equal
from numpy.testing import assert_array_almost_equal as assert_aequal

from impact import Sequence

import os

testdir = os.path.dirname(__file__)

if np.__version__ < '1.12.0':
    calc_prec = 8
else:
    calc_prec = 10


class TestSolenoid1(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_solenoid'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 1.30453956427574,  1.70097634169893, -0.36274921622186,
                         -0.95661373251567, -0.9585393113208 ,  1.32233322144989,
                          1.1655993738424 ,  0.79401111236096,  0.94237284641501,
                         -1.87077312659179,  0.67456890375788,  1.45893510816629]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-1.35895700833264,  0.06910942639346, -3.73107297718927,
                         -3.51931938788647,  3.20471534580495,  0.06851820477388,
                          1.96202863078969,  0.82563918877382, -1.02431346027454,
                         -1.27319733725411,  1.65048729139217, -1.14870689447157]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([-0.75921474578969,  0.07878244391005,  1.22377390297446,
                         -0.87508653029216, -1.67915461947436,  0.9548287825255 ,
                         -1.55213829348015,  1.50901368436332, -2.15674860338446,
                         0.06889642499056,  0.08808255812286,  2.34477928553541]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 3.51384019839184,  2.55155511878505,  0.71033164405217,
                         -2.48472216899002, -1.27116725211923,  0.14681495258728,
                         -3.12342243912057,  3.45319845074631, -3.22429734537948,
                         -1.68811250233734, -1.84158233686758, -0.15482939049823]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([13.06374646465095,  11.09622710173412, -11.28226526138461,
                        -11.33542148314469,  -1.89554838402001,  -9.18177970187455,
                          1.01767325205487,  -2.44705153882094,  -7.07748398778015,
                         13.85311347898182,   7.44048646781677,  -1.60358526033198]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hx3rd('mm')[-1], 0.12663527470631925,
            decimal=calc_prec)

        assert_equal(self.sq.hxp3rd('mrad')[-1], 0.7019179052655159,
            decimal=calc_prec)

        assert_equal(self.sq.hx4th('mm')[-1], 1.6498986258707298,
            decimal=calc_prec)

        assert_equal(self.sq.hxp4th('mrad')[-1], 2.7082157564083862,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hy3rd('mm')[-1], 0.30460030240673952,
            decimal=calc_prec)

        assert_equal(self.sq.hyp3rd('mrad')[-1], 0.48426294209855353,
            decimal=calc_prec)

        assert_equal(self.sq.hy4th('mm')[-1], 1.6607134345690129,
            decimal=calc_prec)

        assert_equal(self.sq.hyp4th('mrad')[-1], 2.7402216011331451,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hrefphi('rad')[-1], 36.059358888718478,
            decimal=calc_prec)

        assert_equal(self.sq.hz3rd('deg')[-1], 7.8075568968200306,
            decimal=calc_prec)

        assert_equal(self.sq.hzp3rd('MeV')[-1], 0.000930066136354909,
            decimal=calc_prec)

        assert_equal(self.sq.hz4th('deg')[-1], 8.5822403470055217,
            decimal=calc_prec)

        assert_equal(self.sq.hzp4th('MeV')[-1], 0.0034720297053829527,
            decimal=calc_prec)


class TestSolenoid2(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_solenoid'))
        self.sq.Flaginteg = 2
        self.sq.Flagdiag = 2
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 1.29140968221311,  1.69614774683912, -0.35773939668474,
                         -0.94884467122492, -0.95678637375909,  1.31844464445247,
                          1.16086351653277,  0.79207027450412,  0.94174361426912,
                         -1.86932171442096,  0.68020835645062,  1.4564483810339 ]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-1.36913193813676,  0.05302768946672, -3.74220796077512,
                         -3.52763522652609,  3.2179177031199 ,  0.06750952230354,
                          1.95146987165583,  0.82694996609426, -1.02599730218082,
                         -1.24783293053045,  1.62820004445475, -1.15721583685075]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([-0.75119970758804,  0.07823091547886,  1.21918429163786,
                         -0.87032794299198, -1.67357180709469,  0.95242874356605,
                         -1.55307313484519,  1.50328834936184, -2.14734441890544,
                          0.07310156148486,  0.07502620057443,  2.33740873189254]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 3.54703408474391,  2.55805931067091,  0.71038442966082,
                         -2.49504571965798, -1.26910185392616,  0.14639497374149,
                         -3.12362116990907,  3.45818568679786, -3.23288322369333,
                         -1.68755970307049, -1.86045051907427, -0.16292057966406]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.1177088667383 ,  11.12886746174933, -11.20869019813506,
                         -11.26061597444544,  -1.87778521117957,  -9.14106557749437,
                           1.03940836918399,  -2.4070289413435 ,  -7.0142437987936 ,
                          13.9318393451756 ,   7.53695146483905,  -1.57784665369527]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 0.69623021444675215,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 0.79259809032695316,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 0.7104853716987602,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 0.7622109347606528,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 0.89699410125984691,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 0.97501115240722969,
            decimal=calc_prec)


class TestQuadrupole1(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_quadrupole'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 2.64571796995758,  1.94915952366665, -1.06582781426569,
                          0.01424329216466,  0.15678070687138,  0.42335363108736,
                          2.15072567432004, -0.38366679764744,  2.85721753706302,
                         -2.03146417012354,  0.20552247552404, -0.70182898674707]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-3.81217505512319, -2.24951068216865, -1.66035393754512,
                         -0.57815840706285,  2.15382363008224, -0.26346382826417,
                          1.35863260802385, -1.05633337440297, -0.47996962022663,
                          1.26519054679985,  1.61955870100351, -0.12945206665692]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 0.26161802186654,  1.29980032527875,  0.99471152310319,
                         -1.11412450396397, -2.42654010541361,  1.93848365396276,
                         -0.2808913291265 ,  1.66112737508227, -0.76613206773545,
                         -1.29034119203369,  0.66709035215037,  3.39014089192646]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.77562053368578,  3.88615498754573,  0.08580280234303,
                         -5.47612336589203, -3.2920730771998 ,  3.61199537384651,
                         -1.23615399888016,  5.49981563118721, -3.86792245460547,
                         -4.04362174385235,  1.04962439911681,  5.35153475741837]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.06374646465103,  11.09622710173413, -11.28226526138472,
                         -11.33542148314466,  -1.89554838402   ,  -9.18177970187453,
                           1.01767325205486,  -2.44705153882096,  -7.07748398778016,
                          13.85311347898172,   7.44048646781681,  -1.60358526033197]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hx3rd('mm')[-1], 0.4853743791249282,
            decimal=calc_prec)

        assert_equal(self.sq.hxp3rd('mrad')[-1], 0.27902056529063884,
            decimal=calc_prec)

        assert_equal(self.sq.hx4th('mm')[-1], 2.3980956744725588,
            decimal=calc_prec)

        assert_equal(self.sq.hxp4th('mrad')[-1], 2.5101559289602391,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hy3rd('mm')[-1], 0.51802078432323073,
            decimal=calc_prec)

        assert_equal(self.sq.hyp3rd('mrad')[-1], 0.79956731768988532,
            decimal=calc_prec)

        assert_equal(self.sq.hy4th('mm')[-1], 2.0614666637317822,
            decimal=calc_prec)

        assert_equal(self.sq.hyp4th('mrad')[-1], 4.0015621030311195,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hrefphi('rad')[-1], 36.059358888724866,
            decimal=calc_prec)

        assert_equal(self.sq.hz3rd('deg')[-1], 7.8075568968200333,
            decimal=calc_prec)

        assert_equal(self.sq.hzp3rd('MeV')[-1], 0.000930066136354909,
            decimal=calc_prec)

        assert_equal(self.sq.hz4th('deg')[-1], 8.5822403470055217,
            decimal=calc_prec)

        assert_equal(self.sq.hzp4th('MeV')[-1], 0.0034720297053829527,
            decimal=calc_prec)


class TestQuadrupole2(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_quadrupole'))
        self.sq.Flaginteg = 2
        self.sq.Flagdiag = 2
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 2.64343384309621,  1.94877738352795, -1.05624231007445,
                          0.01565305592211,  0.15519069589324,  0.422251640445  ,
                          2.14920241096307, -0.38093917419594,  2.8465410053912 ,
                         -2.03589632722571,  0.21339653987824, -0.69972524472965]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-3.8416252530039 , -2.27064505056413, -1.67230508410299,
                         -0.58079181933164,  2.15714827615586, -0.26386667583942,
                          1.34692389618029, -1.05938389229563, -0.47620890932269,
                          1.29700657568754,  1.6068938298617 , -0.12899363819994]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 0.26576824041402,  1.30325597150658,  0.9950522098926 ,
                         -1.10414585370504, -2.4219838151339 ,  1.93297600677873,
                         -0.280789066918  ,  1.65472550675139, -0.76001170189482,
                         -1.29727075555348,  0.66693999558316,  3.38283576680319]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.78157527252281,  3.88801158555155,  0.08722580631622,
                         -5.47526602031416, -3.28685268825439,  3.61188229524386,
                         -1.23732045721051,  5.50069438373357, -3.86798624523017,
                         -4.04711118541872,  1.04459829828541,  5.34814904338791]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.12201781353213,  11.13425991687263, -11.2076324993919 ,
                         -11.26234753261387,  -1.87343024837573,  -9.13642110285147,
                           1.04638626638895,  -2.40604274958994,  -7.00077710915507,
                          13.93966559511062,   7.53600487217484,  -1.5632892472822 ]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 0.722751477554424,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 0.7778781985625977,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 0.69619024342531466,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 0.74122243110461006,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 0.89716263626010229,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 0.97551261698267044,
            decimal=calc_prec)


class TestMultipole(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_multipole'))
        self.sq.Flaginteg = 2
        self.sq.Flagdiag = 2
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 3.69643969919493,  2.08073381104854,  3.81993771925551,
                          0.06075908579783, -3.21546002985317,  3.83888768363449,
                         -6.04000597126535,  3.65627683399832,  2.39704040934943,
                          2.31059810955214, -1.24593341819737, -0.99019713920465]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([ 1.0858785642867 ,  2.21929058777399,  3.43707235997784,
                          0.01867279215298, -0.42666754147913,  4.12282610383794,
                         -6.87584727792029,  2.57227835604426,  1.67183898581566,
                          2.43229113210694, -0.92355767170876, -0.67800019580131]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 2.59302184085872, -2.29021994821892, -1.57564544571876,
                          1.37624908696276, -0.92777260855917, -1.44080543606342,
                         -2.62712480530959,  0.23190212333123, -2.04198745024786,
                          2.22828888973938, -2.90391698366741, -3.326120121202  ]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.14186173080088, -0.49540397311681, -1.04511176664779,
                          0.89861348319655, -1.28647376556505,  0.79235706730794,
                         -2.48650252896164,  0.07050943011481, -2.0158102524782 ,
                          0.08147167098356, -2.22025407974909, -5.73632957947373]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([-1.43587522317207,   4.20039945418675,   3.42249591562093,
                         -2.41830354837589,  -3.75660373660545,  -8.45116573204095,
                         -6.62254011353651, -12.07127593271161,   6.71716536095402,
                          2.12046508354474, -11.78026552121053,   4.32288395037179]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00012259279207,  0.00172016288607, -0.00043629701427,
                         -0.00225625258461,  0.00041177462831, -0.00418561859952,
                         -0.00063754699175, -0.00229777346431,  0.0023261506228 ,
                          0.0024280660878 , -0.00579278118418, -0.00087176873912]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 1.353991171048026,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 1.858662403850537,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 1.4005179085151367,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 2.3437687982433153,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 1.9138849077566304,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 2.6517175869158405,
            decimal=calc_prec)


class TestHSolenoid(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_hsol'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 1.04828965382912,  0.43814503227163, -0.48529699886845,
                         -0.09704169435953,  0.62531040809459, -0.47230236737247,
                          0.13953589556211, -0.04904059031167,  0.36982163194242,
                         -0.3448714702763 , -0.40091434103957, -1.24151993050793]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-1.28006618124948, -1.61731741189933, -1.57612726897918,
                          2.45689067059164,  2.45209635392452, -1.15154384780415,
                          3.32398933435459, -3.59712969421118,  4.07136054555413,
                          0.98551816564889,  0.95210572148211, -2.41237984772909]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 0.95936183259751,  0.86352989271436,  0.57947680138395,
                          0.25685204056859, -1.11151695365799,  0.62328805049163,
                          0.1931699802868 ,  0.22365705852054,  0.69878836287893,
                         -0.71761473090774,  0.01751109026295,  0.91421989051112]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 0.52280878857584,  1.75910521368662, -2.55684572842211,
                         -2.88188043316047,  0.91841028718955,  1.3834060497781 ,
                          2.46843450879344,  1.15281833815194,  0.59807459685911,
                         -2.65511893758177,  1.68902476705229,  0.78263471015308]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 26.85615417769753,  21.07513214008518, -26.70725647118169,
                         -26.07211629939053,  -6.92504717768875, -21.36729522086478,
                          -0.45648440415998, -11.80139173787173, -20.06451674979765,
                          30.68872121905219,  26.77974367620344,  -9.98855528886285]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 0.35045318610364967,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 0.36061415824469,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 0.34021012814165924,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 0.3808655027812024,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 0.90464056480290078,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 0.96116641812963555,
            decimal=calc_prec)


class TestEQuad(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_equad'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([  6.51593662801959,   9.73333885285597, -26.22508441116042,
                          -3.91943545193621,  18.44193366349377,   3.91755044442947,
                          41.05113385048397, -13.10167404662177,  35.71822257083478,
                         -17.08865623624474,  14.50822852459595, -10.7526952380145 ]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-14.93225711945   , -17.64707675905338,  16.57409472315014,
                           2.39873130164016, -16.84744418772407,  -3.01784605198792,
                         -44.39055715623025,  10.32647238845647, -25.93481072102547,
                          41.64336654499222, -39.51246198615156,   9.1031834112367 ]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ -6.46503592417532, -10.39948566788358,   0.72523628945953,
                         -11.83498079789815,  -3.41404211286136,   7.51566126979727,
                          -0.94839876005073,   9.64389884991116,  -7.8963015677347 ,
                          27.77569146571359, -12.82628633025729,   8.70051943662307]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([-18.39580220653528, -35.26010463747019,  -2.01151441347161,
                          -8.63090202074068,   9.02347233160975,   0.84593343436286,
                           2.0373188312178 ,   1.54534710940918,  -5.05662807355479,
                          69.50043795091258, -29.98321562697464,  -5.96966930483419]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 5.37982361150579e+01,   1.04054036688746e+02,
                         -5.24157777116249e+01,   6.25026293875556e+01,
                         -1.24376510886170e+02,  -1.75681593279754e+02,
                          9.68554135236693e-02,   2.51519106042038e+01,
                          3.41729048904977e+01,  -7.02711042671624e+01,
                         -3.07092098572332e+01,  -1.29203421943517e+02]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 30.992275857067796,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 48.470838832613339,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 27.664716018091042,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 57.658815921343624,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 25.881748597206148,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 30.41687621878517,
            decimal=calc_prec)


class TestDipole(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_dipole'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([-23.22215241487229, -15.38564458540654,  30.56882920605662,
                          29.71800280499356,  14.29835728435402,  25.17350554248746,
                           6.97935563394763,  20.12996948554333,  25.67196361648909,
                         -23.3498301926287 , -28.61355075150419,  18.86700382112297]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-144.87632241376505, -102.60610668939259,  224.88118827147287,
                          209.43490141429464,   92.68566732475392,  171.66077541902877,
                           25.50721886074303,  141.55180706932418,  155.7048995755728 ,
                         -106.79487265365589, -148.89721452271039,  134.84263167173521]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 1.30660413282733,  3.3049791764534 ,  1.45539564602175,
                         -4.54348247391579, -4.46859713481191,  4.4419892485426 ,
                         -1.01298550650462,  5.21418565855902, -3.09350005139864,
                         -3.67720549282977,  0.99226146446356,  7.15529373018798]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 2.91475044573253,   4.06050054676938,   1.71900414724068,
                        -25.63911698742061,  -8.38482255841029,  15.51004662716883,
                         -3.88569959847747,  22.75962269467939, -16.73214748393877,
                         -3.61605900505007,  -2.08906697794286,  19.93048036270339]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 44.26352625418364,  155.71212554387185,  105.08622150211926,
                         148.53008804387306,  134.97595261795959, -119.11854331958172,
                         -35.91030563554153,  -38.91272943974511,  -93.42860825171475,
                          96.77194865662064,   -1.53465680972172,   -8.47781431800079]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 43.624448313676311,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 59.135723871038138,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 1.6002390427481115,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 1.8757594033845697,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 26.385757438552069,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 33.17707869506782,
            decimal=calc_prec)


class TestEDipole(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_edipole'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([  1.94981788506178,   5.23778375100127,  -5.94759095843387,
                           9.72253376854368,  18.70539607552703,  13.57520355804367,
                          34.63855650070409,  -0.97575023582377,  37.94301701375338,
                         -29.38945003991655,   2.64595034128239,  -0.29332168586139]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([ 41.54048817899033,  31.16463801111049, -93.17427083591184,
                         -78.17228233766429, -21.31347282295064, -61.32330234058725,
                           4.65532007466343, -56.78666354173257, -49.06757274928654,
                          35.49371804006238,  58.6229778697524 , -51.00819633084841]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 28.41964516969838,  62.71569537564228,   4.68587478301732,
                         -67.66007603052945, -51.03608544789021,  50.09374651964571,
                         -17.18967378206988,  73.005376551632  , -48.35868671431976,
                         -70.75803355563278,  20.68287806260193,  78.52046170991443]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.58875886651737,  3.54296630934073,  0.11732958967841,
                         -4.97625074852837, -3.08330016788748,  3.33561049581566,
                         -1.12239748647181,  5.01973928221869, -3.51729926226983,
                         -3.67509682606576,  0.98224309554126,  4.97319892733316]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([  -8.03456206160189,   91.2138052094256 ,   45.57313926423456,
                           72.01628281929924,   35.10167091292111,  -99.57421992110602,
                            3.44540309230721,  -76.62950376139089, -113.1576698669498 ,
                          -58.89786850380105,   -0.63689973504915,  140.26287182952888]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 38.730857655571761,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 44.097493011725781,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 1.4098923429798638,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 1.8070572922277492,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 26.892991180925677,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 30.957709174083128,
            decimal=calc_prec)


class TestCCL(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.subdir = testdir
        self.sq.read(os.path.join(testdir, 'lattice_sample_ccl'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 2.09124007815674,  1.72846776130466, -1.89864746897705,
                         -0.1703900252042 ,  0.90993860087755,  0.45332679545995,
                          3.189501045709  , -0.83562705508391,  3.47748783983558,
                         -2.16598538737413,  0.80463998592595, -0.93486221955825]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([ 1.58092801900652,  1.66017841978463, -3.45225096458677,
                         -0.47675478093592,  2.18038781380123,  0.56710192796065,
                          5.2285794434673 , -1.64564217646758,  4.95157750161682,
                         -2.67955595999895,  1.8051116093382 , -1.4305406337309 ]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 0.64923639193979,  2.03581098348902,  0.88346598903661,
                         -2.24024420141153, -2.86359068790756,  2.52307871739754,
                         -0.53348376763407,  2.71957663933645, -1.56547970402082,
                         -2.06974272514271,  0.82644683788687,  4.16759427029099]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.33289053849591,  3.19862719448808,  0.48145117189634,
                         -4.22348956989205, -3.26611173553738,  3.29707642539356,
                         -0.96707626659598,  4.46828781979955, -2.97718366573985,
                         -3.30338693334726,  0.99908316925916,  5.06878518098592]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([  9.69460034822973,   7.67887472396391, -13.23992752040479,
                         -13.1984125450551 ,  -4.59180218231297, -11.21372383825554,
                          -1.8238492565765 ,  -5.47234187052865,  -9.60769555170925,
                          10.68599179396361,   5.66982020820139,  -4.70398897999305]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00044706255706, -0.00037085873206, -0.00411484989013,
                         -0.00388241642374, -0.0026438133429 , -0.00347711578799,
                         -0.00208612069062, -0.00393917830151, -0.00420722366784,
                          0.00124422133087,  0.00324988425908, -0.00380200072676]),
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hrefphi('deg')[-1], 2288.9482722490438,
            decimal=calc_prec)

        assert_equal(self.sq.hrefeng('MeV')[-1], 0.50861113995656437,
            decimal=calc_prec)


class TestRFQ(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.subdir = testdir
        self.sq.read(os.path.join(testdir, 'lattice_sample_rfq'))
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([-1.01130706534412, -0.37830769605   , -1.53285919403421,
                         -0.31951820509195,  1.36976290426189,  0.08187083080399,
                          2.00009970702803, -0.85709334307343,  1.28495325676515,
                         -0.30071456725335,  1.04931170823242, -0.44392360766699]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-5.97120737464802, -3.76390977003614, -2.40957986181084,
                         -0.9240983935969 ,  2.99555517170163, -0.51909587231069,
                          0.99966155464744, -1.38056056976332, -1.54981660842168,
                          2.65343340559003,  2.02122667449886,  0.16520722563769]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([-0.23791849774732,  1.11738738544053,  1.70264227302648,
                         -0.3084481895762 , -3.50427787834707,  2.3171655138206 ,
                         -0.08850204893371,  1.20153020491887, -0.10708659172731,
                         -1.03546500795023,  0.91562895296627,  4.74341726847562]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 0.37910064744545,  3.63537466599411,  1.50680432050787,
                         -3.43649385221196, -6.08867674685399,  4.22057092274642,
                         -0.73990213228862,  4.32826131947718, -2.19021209480344,
                         -3.60593916877981,  1.81361970072158,  9.04476430084141]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 12.91320670477946,  11.20153931294174, -10.92391682069154,
                         -10.93259261482697,  -1.86310553439684,  -8.9837635403031 ,
                           1.34288159141171,  -1.82066983798853,  -6.38893001107828,
                          13.75062570443962,   6.52053306575984,  -1.33362579606821]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00431235628048,  0.00330527469154, -0.00492254409138,
                         -0.00464594642525, -0.00180926466932, -0.00396737962473,
                         -0.00016519607858, -0.0025758029467 , -0.00371974381772,
                          0.0052851142    ,  0.00537363651958, -0.00286420081117]),
            decimal=calc_prec)

    def test_xhist(self):
        assert_equal(self.sq.hxepsnp('mm-mrad')[-1], 0.7049721423348556,
            decimal=calc_prec)

        assert_equal(self.sq.hxepsnf('mm-mrad')[-1], 0.7290626412264724,
            decimal=calc_prec)

    def test_yhist(self):
        assert_equal(self.sq.hyepsnp('mm-mrad')[-1], 0.8215439762440437,
            decimal=calc_prec)

        assert_equal(self.sq.hyepsnf('mm-mrad')[-1], 1.007217636360946,
            decimal=calc_prec)

    def test_zhist(self):
        assert_equal(self.sq.hzepsnp('mm-mrad')[-1], 0.91182066494684155,
            decimal=calc_prec)

        assert_equal(self.sq.hzepsnf('mm-mrad')[-1], 0.97220731811963967,
            decimal=calc_prec)


class TestErrors1(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.subdir = testdir
        self.sq.read(os.path.join(testdir, 'lattice_sample_errors'))
        self.sq.Flaginteg = 1
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([-0.88190464384265, -0.05262109864056, -1.5387981131768 ,
                         -1.53823814604199,  2.65545405384276,  0.1436913320707 ,
                          0.97944271324506,  1.00577665406734, -0.90765654364653,
                          0.4001759478918 ,  1.23303305928694, -0.50791190029346]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-2.31573680576508, -0.74439224258695, -5.49800388038717,
                         -3.00062204238222,  3.22093068516848, -0.38490131521419,
                          4.25549032412704, -1.8806367234395 ,  1.69365411525029,
                         -2.41181889025232,  2.09439532441964, -2.09941901420805]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 5.37488169168164,  3.27806882575195,  1.45565720891173,
                         -0.01123953180852,  1.1669482765012 ,  0.02211528656738,
                         -1.53105208221857,  3.22334907650993, -0.65909598210745,
                          0.05977577711299, -1.36300242404338, -1.41582398677685]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 17.99401845966301,  11.3922271561097 ,   3.38538757452654,
                          -3.82801989712956,   1.96338674639036,  -0.27452786304184,
                          -7.41481542867255,  12.05363495942855,  -5.59158634475052,
                          -2.61708150463341,  -5.89778413712969,  -4.91445207594165]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.06374646465101,  11.09622710173412, -11.28226526138469,
                         -11.33542148314467,  -1.89554838402001,  -9.18177970187453,
                           1.01767325205486,  -2.44705153882095,  -7.07748398778016,
                          13.8531134789817 ,   7.44048646781681,  -1.60358526033198]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)


class TestErrors2(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.subdir = testdir
        self.sq.read(os.path.join(testdir, 'lattice_sample_errors'))
        self.sq.Flaginteg = 2
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 2.51027910588845,  1.1392628844616 , -0.72802577599119,
                          0.66687474826932,  3.26428147746699, -0.78168151013661,
                          1.17026954195615,  0.36261875596219,  1.43512791305321,
                          0.67211264912859, -0.06529355377023, -2.83167324553312]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([ 2.90118205583724, -0.17584128876561, -2.48492933195732,
                          3.36010137086321,  7.17706715151451, -3.25383916282909,
                          3.58402005758919, -2.98546558714586,  5.16485683830508,
                          2.12528270872834, -0.35783487939354, -7.88483495906809]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 4.66318940126431,  3.46006201548459,  1.30649414271879,
                          0.54487063233113, -3.8876585282619 ,  1.43723187160898,
                          0.44372116466278,  0.44908147272066,  2.76487698326069,
                         -3.43209465962603, -0.94388036226442,  1.60761005280535]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 16.47651364021324,  12.39083511765802,   2.09952739539513,
                          -0.38499205741948, -15.71589389446372,   4.52782824731897,
                           1.87771092263889,   0.21870413070678,  10.09957573777416,
                         -15.76537532029055,  -4.18222264692257,   4.70571653941911]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.27152871705839,  11.18698685475727, -11.19324346375127,
                         -11.2332733383011 ,  -1.71844252825263,  -9.1273684686828 ,
                           1.0563466469243 ,  -2.38889191268388,  -6.93734200991122,
                          14.02656672274776,   7.55135009946182,  -1.50329168520634]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)


class TestAtoBrun(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_quadrupole'))
        self.sq.Flaginteg = 2
        self.sq.Flagdiag = 2
        sample_dist = np.asarray([
                  [ 0.00087578507661,  0.00091078332293, -0.0018712159038 ,
                   -0.00026020799587,  0.00118246379637,  0.00030477351588,
                    0.00282965108647, -0.00089148165433,  0.00267364758424,
                   -0.00146200173653,  0.00097814405549, -0.00077424907954],
                  [-0.00494130732326, -0.00345004031815,  0.00090239308056,
                   -0.00021646690888,  0.00051744243845, -0.00067269065739,
                   -0.0024990996337 ,  0.00016867110966, -0.00409799768355,
                    0.00322823520659,  0.00026841383209,  0.00092238844105],
                  [ 0.00068673950816,  0.00200932643696,  0.00071693770767,
                   -0.00233188546719, -0.00262497906649,  0.00237609108078,
                   -0.00054685401102,  0.0027194616461 , -0.0016325858362 ,
                   -0.00204748915982,  0.00076273234598,  0.00388112588477],
                  [  2.37815957960549e-03,   5.07795882780603e-03,
                    -5.64848361793523e-05,  -7.28076063642505e-03,
                    -4.07035024099685e-03,   4.59224412285567e-03,
                    -1.63839395253902e-03,   7.21763159965857e-03,
                    -5.14621336651991e-03,  -5.29456063342144e-03,
                     1.31636604043967e-03,   6.71344928482563e-03],
                  [ 16.47437863816542,  13.57115798958197, -15.18582917164655,
                   -15.07109848775285,  -3.16198632877674, -12.25633166125457,
                     0.63424334663149,  -4.81132093322841, -10.37192773644588,
                    18.01084657422632,  12.2082435143298 ,  -3.718258190197  ],
                  [ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                   -0.00422235815026, -0.00143142851249, -0.00347510223932,
                   -0.00043338286052, -0.00267228460453, -0.00372364136289,
                    0.00469939938589,  0.00538889679139, -0.0023901709033 ],
                  [  1.40000000000000e-10,   1.40000000000000e-10,
                     1.40000000000000e-10,   1.40000000000000e-10,
                     1.40000000000000e-10,   1.40000000000000e-10,
                     1.40000000000000e-10,   1.40000000000000e-10,
                     1.40000000000000e-10,   1.40000000000000e-10,
                     1.40000000000000e-10,   1.40000000000000e-10],
                  [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
                  [ 1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9., 10., 11., 12.]
            ])
        self.sq.distribute(particles=sample_dist, Energy=510000, Mass=0.5*931494320.0)
        self.sq.run(start='line', end='line')

    def test_xdist(self):
        assert_aequal(self.sq.getx('mm')[0:12],
            np.asfarray([-1.87807053975634, -0.99625383538719, -1.45005672954784,
                         -0.39325213717372,  1.52659940347147, -0.05874958646546,
                          1.55099900790456, -0.83699606726688,  0.49651727033443,
                          0.29879168795617,  1.17412293419319, -0.29123741674435]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp('mrad')[0:12],
            np.asfarray([-2.79120161418404, -2.2740322852498 ,  2.35607241123947,
                          0.19769584444614, -1.07749156272624, -0.58317157469549,
                         -3.98949849470368,  1.02334975739775, -4.43233158608403,
                          2.78097994115145, -0.95477227081026,  1.18254230253207]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety('mm')[0:12],
            np.asfarray([ 1.5707093445033 ,  3.79912804538029,  0.56730125591631,
                         -4.98268888168664, -3.8767579095805 ,  3.88970966738479,
                         -1.14100337684206,  5.27828968548368, -3.50962110695688,
                         -3.9244517015006 ,  1.18444717486294,  6.0236457769403 ]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp('mrad')[0:12],
            np.asfarray([ 3.20378088568751,  6.93023473075426,  0.04750698801249,
                         -9.85046704036899, -5.72639021460144,  6.36490237037263,
                         -2.21985626657191,  9.83250777574486, -6.9597197637261 ,
                         -7.21680021769059,  1.83536109739795,  9.37160840123537]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz('deg')[0:12],
            np.asfarray([ 20.41520750895452,  16.43195256855383, -19.61957175757439,
                         -19.29486385564516,  -4.5974570138807 , -15.74406971001769,
                           0.20180592930028,  -7.48199885866944, -14.09750551319147,
                          22.81732879688148,  17.70530128535936,  -6.10480274575396]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp('MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)


class TestLatticeGen(unittest.TestCase):
    def setUp(self):
        self.sq = Sequence()
        self.sq.read(os.path.join(testdir, 'lattice_sample_solenoid'))
        dr = [0.1, 4, 20, 0, 0.02]
        q1 = {'length': 0.1, 'seg': 50, 'step': 20, 'type': 'quadrupole', 'B2':   8.0, 'aper': 0.07}
        q2 = {'length': 0.1, 'seg': 50, 'step': 20, 'type': 'quadrupole', 'B2': -24.0, 'aper': 0.07}
        q3 = {'length': 0.1, 'seg': 50, 'step': 20, 'type': 'quadrupole', 'B2':  24.0, 'aper': 0.07}
        self.sq.construct([dr, q1, dr, q3, dr])
        self.sq.insert(4, q2)
        self.sq.insert(5, dr)
        self.sq.insert(0, [0.0, 0, 100, -2, 0.0, 1])
        self.sq.distribute()
        self.sq.run()

    def test_xdist(self):
        assert_aequal(self.sq.getx(100, 'mm')[0:12],
            np.asfarray([ 2.64571796995758,  1.94915952366665, -1.06582781426569,
                          0.01424329216466,  0.15678070687138,  0.42335363108736,
                          2.15072567432004, -0.38366679764744,  2.85721753706302,
                         -2.03146417012354,  0.20552247552404, -0.70182898674707]),
            decimal=calc_prec)

        assert_aequal(self.sq.getxp(100, 'mrad')[0:12],
            np.asfarray([-3.81217505512319, -2.24951068216865, -1.66035393754512,
                         -0.57815840706285,  2.15382363008224, -0.26346382826417,
                          1.35863260802385, -1.05633337440297, -0.47996962022663,
                          1.26519054679985,  1.61955870100351, -0.12945206665692]),
            decimal=calc_prec)

    def test_ydist(self):
        assert_aequal(self.sq.gety(100, 'mm')[0:12],
            np.asfarray([ 0.26161802186654,  1.29980032527875,  0.99471152310319,
                         -1.11412450396397, -2.42654010541361,  1.93848365396276,
                         -0.2808913291265 ,  1.66112737508227, -0.76613206773545,
                         -1.29034119203369,  0.66709035215037,  3.39014089192646]),
            decimal=calc_prec)

        assert_aequal(self.sq.getyp(100, 'mrad')[0:12],
            np.asfarray([ 1.77562053368578,  3.88615498754573,  0.08580280234303,
                         -5.47612336589203, -3.2920730771998 ,  3.61199537384651,
                         -1.23615399888016,  5.49981563118721, -3.86792245460547,
                         -4.04362174385235,  1.04962439911681,  5.35153475741837]),
            decimal=calc_prec)

    def test_zdist(self):
        assert_aequal(self.sq.getz(100, 'deg')[0:12],
            np.asfarray([ 13.06374646465103,  11.09622710173413, -11.28226526138472,
                         -11.33542148314466,  -1.89554838402   ,  -9.18177970187453,
                           1.01767325205486,  -2.44705153882096,  -7.07748398778016,
                          13.85311347898172,   7.44048646781681,  -1.60358526033197]),
            decimal=calc_prec)

        assert_aequal(self.sq.getzp(100, 'MeV')[0:12],
            np.asfarray([ 0.00385496672695,  0.00279736299277, -0.0044121172337 ,
                         -0.00422235815026, -0.00143142851249, -0.00347510223932,
                         -0.00043338286052, -0.00267228460453, -0.00372364136289,
                          0.00469939938589,  0.00538889679139, -0.0023901709033 ]),
            decimal=calc_prec)

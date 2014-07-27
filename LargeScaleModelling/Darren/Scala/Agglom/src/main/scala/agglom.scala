import scala.swing._
import scala.swing.event._
import java.awt.{Graphics2D,Color,BasicStroke}
import java.awt.image.{BufferedImage,WritableRaster}
import scala.math._
import scala.util.Random
import scala.annotation.tailrec

object MySwingApp extends SimpleSwingApplication {

 def top = new MainFrame {
  title="Swing app"
  val button=new Button {
   text="Generate"
  }
  val panel=ImagePanel(500,500)
  contents=new BoxPanel(Orientation.Vertical) {
   contents+=button
   contents+=panel
   border=Swing.EmptyBorder(30,30,10,30)
  }
  listenTo(button)
  var nClicks=0
  reactions+={
   case ButtonClicked(b) =>
    panel.repaint()
  }
 }


}


object ImagePanel {
 def apply(x:Int,y:Int)={
  var ip=new ImagePanel()
  ip.preferredSize=new Dimension(x,y)
  ip.bi=new BufferedImage(x,y,BufferedImage.TYPE_BYTE_GRAY)
  ip
 }
}

class ImagePanel extends Panel {

 private var bi:BufferedImage=null

 override def paintComponent(g: Graphics2D) = {
  g.clearRect(0,0,size.width,size.height)
  //g.setColor(Color.blue)
  // g.plotPoint(400,300)
  val r=Random
  var wr=bi.getRaster()
  var big=bi.createGraphics()
  big.setColor(new Colour("white"))
  big.fillRect(0,0,size.width,size.height)
  agglom(size.width/2,0,100,r,wr)
  g.drawImage(bi,0,0,null)
 }

 @tailrec final def agglom(x0: Int,y0: Int,n: Int,r: Random,wr:WritableRaster): Unit = {
  def adjacent(x: Int,y: Int): Boolean = {
    val byte=wr.getSample(x,y,0)
    true
  }
  if (n>0) {
    @tailrec def wander(x: Int,y: Int): (Int,Int) = {
      val u=r.nextDouble
      val xp=if (u<0.5) {x-1} else {x+1}
      val yp=if ((u>0.25)&(u<0.75)) {y-1} else {y+1}
      val xn=min(max(x,0),size.width)
      val yn=min(max(y,0),size.height)
      if (adjacent(xn,yn)) {
        wr.setSample(xn,yn,0,100)
        (xn,yn)
      } else {
        wander(xn,yn)
      }
    }
    wander(x0,y0)
    agglom(x0,y0,n-1,r,wr)
  }
 }



}




